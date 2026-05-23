"""
ORTHOLOG DISCOVERY PIPELINE

This script:
1. Reads FASTA files using pure Python.
2. Extracts the human target sequence by matching 'PM34_HUMAN' in the header.
3. Runs makeblastdb and blastp using subprocess.
4. Parses BLAST tabular output manually.
5. Collects ortholog sequences into a new FASTA file.
6. Runs MUSCLE via subprocess to generate an alignment.
7. Builds a simple phylogenetic tree using pure Python (UPGMA).
8. Outputs the tree in Newick format.
9. Produces PNG plots using matplotlib + SciPy dendrogram.


"""

import os
import re
import subprocess
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

# =========================================================
# NORMALISATION FUNCTION (FINAL FIX)
# =========================================================

import re

def normalize_header(h):
    """
    Normalise FASTA headers so they match across BLAST, MUSCLE, and alignment.
    Removes MUSCLE suffixes like /1-300 or _1.
    Keeps only the first UniProt-style token.
    """
    token = h.split()[0].strip()
    token = re.sub(r"/\d+-\d+$", "", token)   # remove /1-300
    token = re.sub(r"_\d+$", "", token)       # remove _1
    return token

# =========================================================
# HELPER FUNCTIONS
# =========================================================

def ensure_outputs_folder():
    if not os.path.exists("outputs"):
        os.makedirs("outputs")
        print("[INFO] Created 'outputs' folder.")
    else:
        print("[INFO] 'outputs' folder already exists.")

def check_file(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"ERROR: File not found → {path}")
    if not os.access(path, os.R_OK):
        raise PermissionError(f"ERROR: Cannot read file → {path}")
    print(f"[OK] File is present and readable: {path}")

def run_command(cmd):
    print(f"\n[RUNNING COMMAND]\n{cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print("[ERROR] Command failed.")
        print(result.stderr)
        raise RuntimeError(f"Command failed: {cmd}")
    if result.stdout.strip():
        print(result.stdout)
    return result.stdout

# =========================================================
# FASTA PARSING
# =========================================================

def parse_fasta(filepath):
    with open(filepath) as f:
        header = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    yield (header, "".join(seq_lines))
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
        if header:
            yield (header, "".join(seq_lines))

def write_fasta(records, filepath):
    with open(filepath, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

# =========================================================
# STEP 1 — EXTRACT HUMAN TARGET
# =========================================================

def extract_human_target(human_fasta, target_id_substring="PM34_HUMAN"):
    print("\n[STEP 1] Extracting human target sequence (PM34_HUMAN)...")
    for header, seq in parse_fasta(human_fasta):
        if target_id_substring in header:
            write_fasta([(header, seq)], "target_query.fasta")
            print(f"[OK] Found target sequence: {header}")
            return (header, seq)
    raise ValueError(f"ERROR: No sequence containing '{target_id_substring}' found.")

# =========================================================
# STEP 2 — BLASTP
# =========================================================

def run_blast_for_mammals(mammal_fastas):
    print("\n[STEP 2] Running BLASTP for each mammal FASTA...")
    hits = []
    for fasta in mammal_fastas:
        print(f"\n[INFO] Processing: {fasta}")
        run_command(f"makeblastdb -in {fasta} -dbtype prot")
        out = f"{fasta}_blast.txt"
        run_command(
            f"blastp -query target_query.fasta -db {fasta} "
            f"-outfmt 6 -max_target_seqs 1 -evalue 1e-5 -out {out}"
        )
        if not os.path.exists(out) or os.path.getsize(out) == 0:
            print("[WARN] No BLAST hits.")
            continue
        with open(out) as f:
            line = f.readline().strip()
            if not line:
                continue
            best_hit = line.split("\t")[1]
            hits.append((fasta, best_hit))
            print(f"[OK] Best hit: {best_hit}")
    return hits

# =========================================================
# STEP 3 — COLLECT ORTHOLOGS
# =========================================================

def collect_ortholog_sequences(human_record, hits, output_fasta="orthologs.fasta"):
    print("\n[STEP 3] Collecting ortholog sequences...")
    # Store all records with normalized headers for downstream consistency
    norm_h, norm_seq = normalize_header(human_record[0]), human_record[1]
    all_records = [(norm_h, norm_seq)]
    for fasta, hit_id in hits:
        print(f"[INFO] Searching for {hit_id} in {fasta}")
        target_norm = normalize_header(hit_id)
        found = False
        for header, seq in parse_fasta(fasta):
            if normalize_header(header) == target_norm:
                all_records.append((normalize_header(header), seq))
                print(f"[OK] Added {header}")
                found = True
                break
        if not found:
            print(f"[WARN] Could not find {hit_id} in {fasta}")
    write_fasta(all_records, output_fasta)
    print(f"[OK] Saved to '{output_fasta}'")

# =========================================================
# STEP 4 — MUSCLE ALIGNMENT
# =========================================================

def run_muscle_alignment(input_fasta="orthologs.fasta", output_fasta="orthologs_aligned.fasta"):
    print("\n[STEP 4] Running MUSCLE...")
    try:
        run_command(f"muscle.exe -align {input_fasta} -output {output_fasta}")
    except:
        run_command(f"muscle.exe -in {input_fasta} -out {output_fasta}")
    print(f"[OK] Alignment saved to '{output_fasta}'")

# =========================================================
# STEP 5 — UPGMA TREE
# =========================================================

class UPGMANode:
    def __init__(self, name, members=None, height=0.0):
        self.name = name
        self.members = members if members else [name]
        self.height = height
        self.left = None
        self.right = None

    def is_leaf(self):
        return self.left is None and self.right is None

    def newick(self):
        if self.is_leaf():
            return self.name
        left = self.left.newick()
        right = self.right.newick()
        ldist = self.height - self.left.height
        rdist = self.height - self.right.height
        return f"({left}:{ldist:.4f},{right}:{rdist:.4f})"

def compute_distance_matrix(records):
    """
    Return (labels, distance_matrix) where:
      - labels: list of normalised sequence IDs
      - distance_matrix: n x n numpy array of pairwise distances
    """
    labels = [h for h, s in records]
    seqs = [s for h, s in records]
    n = len(seqs)
    mat = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            mismatches = sum(a != b for a, b in zip(seqs[i], seqs[j]))
            d = mismatches / len(seqs[i])
            mat[i, j] = d
            mat[j, i] = d
    return labels, mat


def upgma(dist):
    clusters = {name: UPGMANode(name) for name in dist}
    while len(clusters) > 1:
        min_d = float("inf")
        pair = (None, None)
        keys = list(clusters.keys())
        for i in range(len(keys)):
            for j in range(i+1, len(keys)):
                a, b = keys[i], keys[j]
                if dist[a][b] < min_d:
                    min_d = dist[a][b]
                    pair = (a, b)
        a, b = pair
        new_name = f"{a}_{b}"
        new_height = min_d / 2
        new_node = UPGMANode(new_name,
                             members=clusters[a].members + clusters[b].members,
                             height=new_height)
        new_node.left = clusters[a]
        new_node.right = clusters[b]
        clusters[new_name] = new_node
        for k in list(clusters.keys()):
            if k in (a, b, new_name):
                continue
            d = (
                dist[a][k] * len(clusters[a].members)
                + dist[b][k] * len(clusters[b].members)
            ) / (len(clusters[a].members) + len(clusters[b].members))
            dist[new_name][k] = d
            dist[k][new_name] = d
        del clusters[a]
        del clusters[b]
        dist.pop(a)
        dist.pop(b)
        for v in dist.values():
            v.pop(a, None)
            v.pop(b, None)
    return list(clusters.values())[0]

# =========================================================
# OPTIMISED UPGMA PLOTTING
# =========================================================

def plot_upgma_tree(linkage_matrix, labels, output_path):
    species_colors = {
        "HUMAN": "darkred",
        "MOUSE": "steelblue",
        "GORGO": "forestgreen",
        "PANTR": "orange",
        "CANLF": "purple"
    }

    plt.figure(figsize=(12, 8))
    dendro = dendrogram(
        linkage_matrix,
        labels=labels,
        orientation='right',
        leaf_font_size=10,
        leaf_rotation=0,
        color_threshold=0.7 * max(linkage_matrix[:, 2]),
        above_threshold_color="gray"
    )

    ax = plt.gca()
    for lbl in ax.get_ymajorticklabels():
        text = lbl.get_text()
        for key, color in species_colors.items():
            if key in text:
                lbl.set_color(color)
                lbl.set_fontweight("bold")

    plt.title("UPGMA Phylogenetic Tree of PMP34 Orthologs", fontsize=14, fontweight="bold")
    plt.xlabel("Evolutionary Distance", fontsize=12)

    max_dist = max(linkage_matrix[:, 2])
    plt.hlines(y=-1, xmin=0, xmax=max_dist * 0.2, colors="black", linewidth=2)
    plt.text(max_dist * 0.1, -1.5, "Scale: 0.2 distance", ha="center", fontsize=10)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

    print(f"[OK] Saved optimised UPGMA tree to '{output_path}'")

# =========================================================
# BUILD TREE
# =========================================================

def build_phylogenetic_tree(aligned_fasta="orthologs_aligned.fasta", out_prefix="outputs/ortholog_tree"):
    print("\n[STEP 5] Building phylogenetic tree (UPGMA)...")

    # Normalise MUSCLE headers once
    records = [(normalize_header(h), seq) for h, seq in parse_fasta(aligned_fasta)]

    # Get labels and distance matrix directly
    labels, mat = compute_distance_matrix(records)

    # Build UPGMA tree using labels + matrix
    # Convert matrix to dict-of-dicts for the existing upgma() function
    dist = defaultdict(dict)
    n = len(labels)
    for i in range(n):
        for j in range(n):
            dist[labels[i]][labels[j]] = mat[i, j]

    root = upgma(dist)
    newick = root.newick() + ";\n"

    newick_path = f"{out_prefix}_UPGMA.newick"
    with open(newick_path, "w") as f:
        f.write(newick)
    print(f"[OK] Saved UPGMA tree to '{newick_path}'")

    # Use SciPy linkage on the same matrix
    condensed = squareform(mat)
    linkage_matrix = linkage(condensed, method='average')

    plot_upgma_tree(linkage_matrix, labels, f"{out_prefix}_UPGMA_pretty.png")


# =========================================================
# MAIN
# =========================================================

def main():
    print("\n======================================")
    print("   ORTHOLOG DISCOVERY PIPELINE START  ")
    print("======================================")

    ensure_outputs_folder()

    human_fasta = "uniprotkb_proteome_UP000005640_9606.fasta"
    mammal_fastas = ["mammal_1.fasta", "mammal_2.fasta", "mammal_3.fasta", "mammal_4.fasta"]

    check_file(human_fasta)
    for f in mammal_fastas:
        check_file(f)

    human_record = extract_human_target(human_fasta)
    hits = run_blast_for_mammals(mammal_fastas)
    collect_ortholog_sequences(human_record, hits)
    run_muscle_alignment()
    build_phylogenetic_tree()

    print("\n======================================")
    print("   ORTHOLOG DISCOVERY PIPELINE DONE   ")
    print("======================================")

if __name__ == "__main__":
    main()

