# Mfuzz-analysis-function-repo.R

## count2tpm 
Convert raw counts into a gene matrix normalized to TPM.
- Input the gene matrix `gene_matrix` and the effective gene length file `gene_efflen`
- Output the gene matrix normalized to TPM.
## MDataLoader
load metabolomics data
- Input the file location `files` and `sheetnames`
- Output the standard preprocessed metabolomic data.
## TDataLoader
Load transcriptomics data.
- Input the file location `expression_path`、`sheetnames`、address to sample information `sample_info_path`、address to gene effective length files `gene_efflen_path`.
- Output the transcriptomic data normalized to TPM.
## KW_test 
perform KW test 
- Input the parameter list `paras`、Input data adress `files`、`sheetnames`、output file name `output`.
- Output the results of the KW test and store them in the corresponding Result file.

## ENSEMBL2SYMBOL
Convert `ENSEMBL id` to `Gene SYMBOL`
- Input the gene expression matrix `gene_matirx` and the reference genome `OrgDb`.
- Output `Gene SYMBOL` data.

## Mfuzz_KW
Mfuzz analysis and plot Mfuzz clustering diagrams.
- Input parameter list`paras`、Input file address `files`、`sheetnames`、output file name `output_name`、p-value threshold `pvalue_threshold`、time-series labels `labels`、The number of Mfuzz cluster `Mfuzz_number`.
- Output Mfuzz clustering results and Mfuzz clustering diagrams.
## mfuzz.plot2
Customized Mfuzz plotting function, adapted from the official mfuzz.plot2 fucntion.