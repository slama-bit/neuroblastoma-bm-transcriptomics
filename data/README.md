# Data

The raw data are not tracked in this repository.

| Item | Value |
|---|---|
| GEO series | [GSE25624](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25624) |
| Title | The bone marrow microenvironment in patients with neuroblastoma: presence of IFN signature and significant down-regulation of CXCL12 expression |
| Platform | GPL4133 (Agilent-014850 Whole Human Genome Microarray 4x44K) |
| Samples | 40 bone marrow samples (GSM629613–GSM629652) |

## Groups

| Code | Group | n | GEO annotation |
|---|---|---|---|
| **I** | Metastatic NB (infiltrated BM) | 9 | patient disease status: metastatic Neuroblastoma, tissue: resident bone marrow cells |
| **NI** | Localized NB (non-infiltrated BM) | 23 | patient disease status: localized Neuroblastoma, tissue: resident bone marrow cells |
| **H** | Healthy young donors | 8 | patient disease status: healthy young donor, tissue: bone marrow |

## Download

To download the series in R:

```r
library(GEOquery)
gset <- getGEO("GSE25624", GSEMatrix = TRUE)[[1]]
saveRDS(gset, "data/GSE25624_gset.rds")
```
