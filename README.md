# NCBI RNA-seq Data Retrieval and Processing

An R-based workflow for retrieving publicly available bulk RNA-seq datasets from NCBI GEO, extracting sequencing metadata, and preparing raw gene count matrices for downstream differential gene expression (DGE) analysis using `edgeR`.

## Overview

The workflow provides functions to:

- Retrieve and inspect GEO Series (GSE) metadata and supplementary files.
- Download NCBI-generated raw gene count matrices (GRCh38.p13).
- Retrieve SRA run metadata (`SraRunInfo.csv`) and sample attributes (`SraRunTable.csv`) using NCBI E-utilities.
- Resolve GEO (GSE) accessions to their corresponding SRA Study (SRP) identifiers.
- Integrate gene counts, sample metadata, and gene annotations into `edgeR::DGEList` objects (main output).
- Validate count matrices and sample correspondence before exporting processed objects.

## Requirements

The workflow is implemented in **R** using the following packages:

- `Biobase`
- `edgeR`
- `GEOquery`
- `rentrez`
- `recount3`
- `purrr`
- `tidyverse`
- `xml2`

## Input Data

The pipeline accepts GEO (`GSE`) accessions for metadata retrieval. NCBI-generated gene count downloads and `DGEList` construction use GSE identifiers.

Example datasets include studies associated with Alzheimer's disease (AD), Huntington's disease (HD), and Parkinson's disease (PD).

## Output

The workflow organizes retrieved data into the following directories:

```text
data_ncbi/
├── annotation/     # Gene annotation file (included)
├── counts/         # NCBI raw gene count matrices
└── metadata/       # SRA run and sample metadata

out_edgeR/
├── dge_raw.rds           # Main `DGEList` object for downstream RNA-seq analysis
└── dge_raw_summary.csv   # Summary metrics of studies contained en the `DGEList` object
```

## Notes

- The workflow currently targets human RNA-seq datasets using the GRCh38.p13 reference annotation.
- Count matrices must correspond to the sample identifiers provided in the retrieved metadata.
- The exported objects contain **raw counts**; low-expression filtering, normalization, and differential expression modeling should be performed separately.
- Some datasets may require additional metadata curation before downstream analysis.

## Author

**Moises A. Rojas**
Laboratory of Genomics and Evolution (MolEvol lab) @CINV 2026
