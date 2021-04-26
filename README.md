# RNA_Seq_Template

This is a template for the downstream analysis of a bulk RNA seq experiment. It requires a count matrix and meta data file and performs a differential expression analysis using DESeq2 and a GSEA enrichment analysis for each contrast defined in the config file. In addition basic QC plots (PCA, coverage, ...) are produced.

The config file needs to be filled in for each run and R packages have to be installed. Container version is in progress.