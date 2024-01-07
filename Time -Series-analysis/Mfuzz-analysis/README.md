# Mfuzz-analysis
Fuzzy c-means clustering, a method used for gene expression data analysis, can be used to identify characteristic expression patterns between different biological samples, grouping them into distinct feature clusters, helping researchers understanding expression changes under different conditions.
## Abstract
`Mfuzz-analysis.R` can be used for analyzing both metabolomic and transcriptomic data modes. It ultimately outputs Mfuzz clustering results and Mfuzz cluster diagrams.

- `Data`：This directory is used to store the raw data related to Mfuzz analysis. The `Transcriptomics` folder is for storing transcriptomic data including the transcriptomic raw data `PDX-read-count.txt`, effective gene length file `gene_efflen.csv`,sample information files `Sample_information.csv`, etc. 
- `Script`：This directory is for storing scripts related to Mfuzz analysis, including the run script `Mfuzz-analysis.R` and the function library `Mfuzz-analysis-function-repo.R`.
- `Result`：The directory is for storing result files, including KW test results, Mfuzz clustering results, and Mfuzz clustering diagrams.