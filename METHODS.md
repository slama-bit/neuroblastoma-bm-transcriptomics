# Methods

## Data acquisition

Gene expression data were retrieved from the NCBI Gene Expression Omnibus (GEO) under accession **GSE25624** using the GEOquery package in R (`getGEO(..., GSEMatrix = TRUE, AnnotGPL = TRUE)`). The annotated series matrix was imported together with GPL-derived probe annotations. Profiling had been performed on the Agilent-014850 Whole Human Genome Microarray 4×44K G4112F (GEO platform **GPL4133**).

Samples were assigned to three groups using the group suffix of the GEO sample titles:

| Group | Description | n |
|---|---|---|
| H_BM | Bone marrow from healthy young donors | 8 |
| NI_BM | Resident bone marrow cells from patients with localized neuroblastoma (no bone marrow infiltration) | 23 |
| I_BM | Resident bone marrow cells from patients with metastatic neuroblastoma (bone marrow infiltration) | 9 |

## Preprocessing and normalization

The distribution of the expression values was checked using quantiles. Values were log2-transformed only if they appeared to be on a linear scale (99th percentile > 100, or range > 50); in that case non-positive values were first set to missing. Between-array quantile normalization was then applied with `normalizeBetweenArrays()` from the limma package.

## Differential expression analysis

Differential expression was tested at the probe level with limma. A linear model with one coefficient per group and no intercept (`~0 + group`) was fitted to the normalized data, and standard errors were moderated with empirical Bayes (`eBayes`). Three contrasts were evaluated:
- metastatic vs healthy (I_BM − H_BM)
- localized vs healthy (NI_BM − H_BM)
- metastatic vs localized (I_BM − NI_BM)

For each contrast, probe-wise log2 fold changes, moderated t-statistics and P values were computed. P values were corrected for multiple testing across all probes with the Benjamini–Hochberg false discovery rate (FDR) procedure.

## Probe-to-gene summarization and DEG definition

Probe-level results were summarized to gene symbols in four steps:
1. Probes with FDR < 0.05 were retained.
2. Probes lacking a gene symbol annotation were removed.
3. For genes represented by more than one probe, the probe with the lowest adjusted P value was kept.
4. The |log2 fold change| ≥ 1 threshold was applied to that representative probe.

Differentially expressed genes (DEGs) were therefore defined as genes with FDR < 0.05 and |log2FC| ≥ 1.

## Visualization

**Heatmaps.** Heatmaps were drawn with pheatmap for the 15 most significant DEGs in the I_BM vs H_BM and NI_BM vs H_BM contrasts. Each heatmap was restricted to the samples of the two groups in that contrast. Expression values were z-scored per gene, and genes and samples were clustered hierarchically (pheatmap defaults: Euclidean distance, complete linkage).

**Volcano plots.** Volcano plots (ggplot2, ggrepel) show the gene-level DEGs of each contrast. The five most significant genes per contrast are labelled.

**Overlap analysis.** The DEG sets were compared with Venn diagrams (VennDiagram). Two core signatures were defined:
- The **disease-associated core** is the set of genes that are DEGs in both I_BM vs H_BM and NI_BM vs H_BM.
- The **metastasis-associated core** is the set of genes that are DEGs in both I_BM vs H_BM and I_BM vs NI_BM.

Genes found in only one contrast were exported as contrast-specific sets.

## Functional enrichment analysis

Each contrast's gene-level DEG list was analysed separately with two tools.

**PANTHER.** Gene Ontology overrepresentation was tested with the PANTHER Overrepresentation Test (PANTHER 19.0, GO-Slim Biological Process, release 2024-06-20). The test used Fisher's exact test with FDR correction, against the default *Homo sapiens* whole-genome reference list. PANTHER was run for the I_BM vs H_BM and NI_BM vs H_BM contrasts.

**g:Profiler.** Functional profiling with g:Profiler (g:GOSt) covered GO (BP, MF, CC), KEGG, Reactome, CORUM and miRNA targets, and was run for all three contrasts. Results with an adjusted P < 0.05 were considered significant. *[Confirm the settings used: organism, multiple-testing method (the default is g:SCS), and whether the query was ordered.]*

## Protein–protein interaction network and hub genes

Protein–protein interaction (PPI) networks were built from the I_BM vs H_BM and NI_BM vs H_BM DEG lists using the STRING database (*Homo sapiens*) and imported into Cytoscape. *[Add the STRING version and the interaction-confidence cut-off, e.g. 0.4.]* Nodes were ranked by degree centrality, and the 30 highest-degree nodes were taken as hub genes. Each hub-gene subnetwork was extracted together with the shortest paths connecting the hubs. Hub status describes network connectivity only and does not by itself show a regulatory or causal role.

## Software and reproducibility

The differential expression, normalization, heatmap, volcano and Venn analyses were run in R using Bioconductor (GEOquery, limma) and CRAN packages (dplyr, ggplot2, ggrepel, VennDiagram, pheatmap). The code is in `scripts/deg_analysis.R`, and package versions are recorded in `results/sessionInfo.txt` when the script is run.

The enrichment analyses (PANTHER, g:Profiler web servers) and the network analysis (STRING, Cytoscape desktop) were done interactively. Their exported results are included in `enrichment/` and `network/`.

Statistical significance was defined as an FDR-adjusted P < 0.05 unless stated otherwise.
