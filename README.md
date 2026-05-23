# Ortholog Discovery and Phylogenetic Analysis of PMP34 (SLC25A17)

This repository contains a complete computational pipeline for identifying orthologous sequences of the peroxisomal membrane protein PMP34 (SLC25A17) across mammalian species and reconstructing their evolutionary relationships using a transparent, dependency minimal Python workflow. The project demonstrates a full comparative genomics workflow suitable for MSc level research in bioinformatics, drug discovery, and molecular evolution.

---

## Abstract

Ortholog identification and phylogenetic reconstruction are fundamental tasks in comparative genomics. This project implements a fully transparent, Python based pipeline for discovering orthologs of the peroxisomal membrane protein PMP34 (SLC25A17) across mammalian species. Using BLASTP for homology search, MUSCLE for multiple sequence alignment, and a custom UPGMA algorithm for tree construction, the workflow quantifies evolutionary divergence and visualises species relationships. The resulting phylogeny reflects established mammalian taxonomy, confirming the evolutionary conservation of PMP34 and demonstrating the utility of lightweight computational pipelines for protein evolution studies.

---

## 1. Scientific Goal

The objective of this project is to investigate the evolutionary conservation of the peroxisomal transporter PMP34 across mammals. PMP34 is a member of the mitochondrial carrier family and plays a key role in peroxisomal metabolite transport. Understanding its conservation provides insight into functional stability, evolutionary divergence, suitability of model organisms, and structural constraints acting on membrane transporters.

This project implements a reproducible computational pipeline to identify orthologs, align their sequences, compute evolutionary distances, and construct a phylogenetic tree.

---

## 2. Workflow Overview

### 2.1 Workflow Diagram

```mermaid
flowchart TD
    A[Human Proteome FASTA] --> B[Extract PMP34 Sequence]
    B --> C[BLASTP Against Mammalian Proteomes]
    C --> D[Collect Top Ortholog Hits]
    D --> E[MUSCLE Multiple Sequence Alignment]
    E --> F[Compute Pairwise Distances]
    F --> G[UPGMA Tree Construction]
    G --> H[Newick Export]
    G --> I[Dendrogram Visualisation]


## 3. Workflow Rationale
The pipeline follows a standard comparative genomics workflow:

Extraction of the human PMP34 sequence from the human proteome FASTA.

Identification of orthologs via BLASTP searches against mammalian proteomes (Camacho et al., 2009).

Consolidation of top scoring orthologs into a single FASTA file.

Multiple sequence alignment using MUSCLE (Edgar, 2004; Edgar, 2024).

Computation of pairwise evolutionary distances.

Construction of a UPGMA phylogenetic tree.

Visualisation of the tree using SciPy’s dendrogram function (Virtanen et al., 2020).

This workflow is intentionally implemented without Biopython to demonstrate algorithmic understanding and ensure full transparency of each computational step.

## 4. Methods
4.1 Data Sources
Protein FASTA files were obtained from UniProtKB for:

Homo sapiens

Mus musculus

Gorilla gorilla gorilla

Pan troglodytes

Canis lupus familiaris

4.2 FASTA Parsing
Custom Python functions were used to read and write FASTA files, ensuring compatibility across proteomes and avoiding external dependencies.

4.3 BLASTP Searches
Each mammalian proteome was converted into a BLAST database using makeblastdb.
The human PMP34 sequence was used as a query in blastp, retrieving the top scoring ortholog per species (Camacho et al., 2009).

4.4 Ortholog Collection
The best BLAST hits were extracted and consolidated into orthologs.fasta.

4.5 Multiple Sequence Alignment
MUSCLE (v5) was executed via subprocess to generate a high quality alignment (Edgar, 2004; Edgar, 2024).
Header normalisation was applied to remove MUSCLE generated suffixes.

4.6 Distance Matrix
Pairwise distances were computed as the proportion of mismatched amino acids between aligned sequences.
This metric reflects evolutionary divergence under the assumption of equal substitution rates.

4.7 UPGMA Tree Construction
A custom UPGMA implementation was used to cluster sequences based on average linkage.
The resulting tree was exported in Newick format.

4.8 Tree Visualisation
SciPy’s dendrogram function was used to generate a publication quality phylogenetic tree with species specific colouring and a distance scale bar (Virtanen et al., 2020).

## 5. Results
5.1 Ortholog Identification
BLASTP successfully identified PMP34 orthologs in all five mammalian species.
Each ortholog exhibited high sequence similarity to the human PMP34, confirming conserved function across taxa.

5.2 Multiple Sequence Alignment
The MUSCLE alignment revealed strong conservation of transmembrane helices and peroxisomal targeting motifs, consistent with PMP34’s essential role in peroxisomal transport.

5.3 Phylogenetic Reconstruction
Figure 1. UPGMA phylogenetic tree of PMP34 orthologs across five mammalian species.  
The tree shows the evolutionary relationships among PMP34 orthologs from human, chimpanzee, gorilla, mouse, and dog. Branch lengths represent the fraction of amino acid mismatches. The primate sequences cluster tightly, while mouse and dog form more distant branches, reflecting earlier divergence.

5.4 Interpretation
Primate Clustering
Human, chimpanzee, and gorilla sequences form a tight cluster with minimal branch lengths, indicating recent divergence and high sequence conservation.

Rodent and Canine Divergence
Mouse and dog sequences form a separate clade with longer branch lengths, reflecting greater evolutionary distance from primates.

Functional Conservation
The overall low divergence values demonstrate that PMP34 is highly conserved across mammals.
This suggests strong purifying selection to maintain peroxisomal transport function.

Biological Implications
PMP34’s conservation supports its essential metabolic role.
Model organisms such as mouse and dog retain sufficient similarity for functional studies.
The phylogeny aligns with established mammalian taxonomy, validating the computational approach.

## 6. Conclusion
This project demonstrates a complete, transparent, and reproducible pipeline for ortholog discovery and phylogenetic analysis using pure Python.
The resulting phylogenetic tree confirms the evolutionary conservation of PMP34 across mammals and illustrates the power of computational genomics for understanding protein evolution.

The workflow is suitable for extension to other protein families, comparative genomics studies, or functional annotation pipelines.

## 7. How to Cite This Repository
If you use this pipeline in academic work, please cite:

Anjali (2026). Ortholog Discovery and Phylogenetic Analysis of PMP34 (SLC25A17). GitHub Repository. Liverpool, United Kingdom.

## 8. References
Edgar, R. C. (2004). MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Research.
Edgar, R. C. (2024). MUSCLE v5: improved alignment accuracy and performance.
Camacho, C. et al. (2009). BLAST+: architecture and applications. BMC Bioinformatics.
Virtanen, P. et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. Nature Methods.