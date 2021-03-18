# This is the analysis script for project #XXXX. This should be sourced before compiling the .Rmd report.

# Load packages
suppressMessages(library("ggplot2"))
suppressMessages(library("tidyr"))
suppressMessages(library("DESeq2"))
suppressMessages(library("pheatmap"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("modules"))

# Set working directory to main.R script location
proj_dir <- dirname(sys.frame(1)$ofile)
# Or select proj_dir for interactive session
#proj_dir <- "~/Desktop/NBIS/Projects/project_5566/"
setwd(proj_dir)

# Load several custom functions
lib <- modules::use("R")

#############
# Load data #
#############

# Set up a report folder, remove existing ones first
res_dir <- paste0(proj_dir, "/Report/")
unlink(res_dir, recursive = TRUE)
dir.create(res_dir, showWarnings = FALSE)

# Load datasets, make sure sample names are the same in data and meta
cts   <- read.csv(paste0(proj_dir, "/data/merged_gene_counts.txt"), sep = "\t", row.names = 1, check.names = FALSE)
cts <- cts[,-1]
rownames(cts) <- gsub(".[0-9]+$", "", rownames(cts))
meta   <- read.csv(paste0(proj_dir, "/data/metadata.csv"), sep = ";", row.names = 1)

# make sure samples are in the same order in metadata and count data
meta <- meta[colnames(cts),, drop=FALSE]
cts    <- cts[,rownames(meta)]
if(!identical(colnames(cts), rownames(meta))){stop("Samplenames do not match between counts and metadata")}

# Make and Save DESeq object
dds <- DESeqDataSetFromMatrix(countData = cts, colData = meta, design = ~ 1)

##############
# Pre Filter #
##############

# Remove all zero genes
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]

###################
# Quality Control #
###################

# Create a QC folder in the Report dir
qc_dir <- paste0(res_dir, "QC/")
dir.create(qc_dir, showWarnings = FALSE)

# Create barplot of the total number of reads per sample
pdf(paste0(qc_dir, "barplot_raw_counts.pdf"))
par(mar=c(7,6,2,2), mgp=c(5,1,0))
barplot(colSums(counts(dds)), 
        col = ifelse(colData(dds)$Status == "Control", "green", "red"),
        ylab = "Number of mapped reads", las =2, cex.names=0.8)
abline(h = mean(colSums(counts(dds))))
dev.off()

# Create boxplot of logged raw counts: are there big differences in the library size?
pdf(paste0(qc_dir, "boxplot_raw_counts.pdf"))
par(mar=c(6,8,2,2))
boxplot(log2(counts(dds)+1), 
        col = ifelse(colData(dds)$Status == "Control", "green", "red"),
        pch = ".",
        horizontal = TRUE,
        las = 1,
        cex.names=0.5,
        xlab = "log2(Counts +1)")
dev.off()

# PCA plot of 100 most variable genes: how is the data structured?
dds   <- estimateSizeFactors(dds)
rld   <- vst(dds)
p1    <- plotPCA(rld, intgroup = "Status", ntop = 100)
ggsave(paste0(qc_dir, "PCA_rlog_top100.pdf"), plot = p1)
p2    <- plotPCA(rld, intgroup = "Status", ntop = 500)
ggsave(paste0(qc_dir, "PCA_rlog_top500.pdf"), plot = p2)
p3    <- plotPCA(rld, intgroup = "Status", ntop = 1000) + geom_label(aes(label = colnames(dds)))
ggsave(paste0(qc_dir, "PCA_rlog_top1000.pdf"), plot = p3)

# Create boxplot of logged raw counts: are there big differences in the library size?
pdf(paste0(qc_dir, "boxplot_normalized_counts.pdf"))
par(mar=c(6,8,2,2))
boxplot(assay(rld), 
        col = ifelse(colData(dds)$Status == "Control", "green", "red"),
        pch = ".",
        horizontal = TRUE,
        las = 1,
        xlab = "Counts (normalized by rlog)")
dev.off()

# Sample Distances: how are the samples related to each other?
sampleDists                 <- dist(t(assay(rld)))
sampleDistMatrix            <- as.matrix(sampleDists)
rownames(sampleDistMatrix)  <- paste(meta$samples, meta$number.of.eggs, sep = "_" )
colnames(sampleDistMatrix)  <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pdf(paste0(qc_dir, "HeatmapDistances.pdf"))
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.off()

##################
# Filter Samples #
##################

# Remove samples if warranted by QC and/or other reasons
dds <- dds[, !colnames(dds) %in% c("F36")]

# PCA plot of 100 most variable genes: did the sample filtering help?
dds   <- estimateSizeFactors(dds)
rld   <- vst(dds)
p1    <- plotPCA(rld, intgroup = "Status", ntop = 100)
ggsave(paste0(qc_dir, "PCA_rlog_top100_Filt.pdf"), plot = p1)
p2    <- plotPCA(rld, intgroup = "Status", ntop = 500)
ggsave(paste0(qc_dir, "PCA_rlog_top500_Filt.pdf"), plot = p2)
p3    <- plotPCA(rld, intgroup = "Status", ntop = 1000) + geom_label(aes(label = colnames(dds)))
ggsave(paste0(qc_dir, "PCA_rlog_top1000_Filt.pdf"), plot = p3)


############
# Heatmaps #
############

# Create a differential expression result folder
sig_dir <- paste0(res_dir, "Signatures/")
dir.create(sig_dir, showWarnings = FALSE)

# Create heatmaps for interferon responses
goList <- gage::readList("data/Human_GOBP_AllPathways_no_GO_iea_March_01_2021_symbol.gmt.txt")
names(goList) <- gsub("%.*", "", names(goList))
interferonSignatures <- goList[grep("INTERFERON", names(goList))]

for(sig in names(interferonSignatures)){
pdf(paste0(sig_dir, sig, "_clustered.pdf"))
lib$plotCustomHeatmap(obj = rld, 
                  plotGenes =  interferonSignatures[sig],
                  anotationColumn = "Status",
                  anotationColor = list(Status = c("Control" = "green", "Fibro" = "purple")),
                  convertToSymbol = TRUE, 
                  convertFromID = "ENSEMBL",
                  db = org.Hs.eg.db, 
                  groupColumns = FALSE)
dev.off()

pdf(paste0(sig_dir, sig, "_grouped.pdf"))
lib$plotCustomHeatmap(obj = rld,
                  plotGenes =  interferonSignatures[sig],
                  anotationColumn = "Status",
                  anotationColor = list(Status = c("Control" = "green", "Fibro" = "purple")),
                  convertToSymbol = TRUE,
                  convertFromID = "ENSEMBL",
                  db = org.Hs.eg.db, 
                  groupColumns = TRUE)
dev.off()
}

###########################
# Differential Expression #
###########################

# Create a differential expression result folder
de_dir <- paste0(res_dir, "Differential_Expression/")
dir.create(de_dir, showWarnings = FALSE)

designList   <- list(Status = list(Model = "~ 0 + Status",
                                   Contrast = list("Fibro_vs_Control" = c("Status", "Fibro", "Control"))))

# Run DE and plotting function over each comparison in the designlist
for (design in names(designList)){
  
  design_dir <- paste0(de_dir, design, "/")
  dir.create(design_dir, showWarnings = FALSE)
  
  # Update design and run DESeq
  dds <- DESeqDataSet(dds, as.formula(designList[[design]][["Model"]]))
  dds <- DESeq(dds)
  
  # Create results for each contrast in design parameter
  for(contrast in names(designList[[design]][["Contrast"]])) {
    
    # Make separate folder for each contrast
    contrast_dir <- paste0(design_dir, contrast, "/")
    dir.create(contrast_dir, showWarnings = FALSE)
    
    # create DE folder in each contrast folder
    diff_exp_dir <- paste0(contrast_dir, "DE/")
    dir.create(diff_exp_dir, showWarnings = FALSE)
    
    # Extract the comparison for plotting purposes
    comp <- designList[[design]][["Contrast"]][[contrast]]

    # Extract results
    res <- results(dds, contrast = comp)

    # Convert ensembl id to symbols
    res$Symbol <- mapIds(org.Hs.eg.db, 
                         keys = rownames(res), 
                         keytype = "ENSEMBL", 
                         column = "SYMBOL")
    res <- res[order(res$padj),]
    write.csv(res, paste0(diff_exp_dir, paste(comp, collapse="_"), ".csv"))
    
    # Create volcano plot
    png(paste0(diff_exp_dir, contrast, "_volcano.png"))
    p1 <- lib$volcano(res, fc = 1, sig = 0.05)
    print(p1)
    dev.off()

    ##################
    # GSEA Enrichment#
    ##################

    # Log2Fc statistic
    ##################
    
    genelist <- sort(setNames(res$log2FoldChange, rownames(res)), decreasing = TRUE)
    
    # genelist names should be entrezid
    names(genelist) <- mapIds(org.Hs.eg.db, 
                              keys = names(genelist), 
                              column = "ENTREZID", 
                              keytype = "ENSEMBL")
    genelist <- genelist[!is.na(names(genelist))]
    
    for(col in c("MF", "BP", "CC")){
    # GO BP Enrichment
    gsea_dir_tmp <- paste0(contrast_dir, "GSEA/Log2FC/GO/", col, "/")
    dir.create(gsea_dir_tmp, showWarnings = FALSE, recursive = TRUE)
  
    lib$visualizeGSEA(geneList = genelist, 
                go_class = col, 
                n_terms = 20, 
                outdir = gsea_dir_tmp, 
                comp=comp,
                collection = "GO",
                org = "hsa",
                orgdb = org.Hs.eg.db,
             msigdb_file = "data/Human_MSigdb_March_01_2021_Entrezgene.gmt")
    }
    
    # KEGG Enrichment
    kegg_dir <- paste0(contrast_dir, "GSEA/Log2FC/KEGG/")
    dir.create(kegg_dir, showWarnings = FALSE, recursive = TRUE)
    
    ib$visualizeGSEA(geneList = genelist, 
             go_class = "", 
             n_terms = 20, 
             outdir = kegg_dir, 
             comp=comp,
             collection = "KEGG",
             org = "hsa",
             orgdb = org.Hs.eg.db,
             msigdb_file = "data/Human_MSigdb_March_01_2021_Entrezgene.gmt")
    
    # MSIGDB Enrichment
    msigdb_dir <- paste0(contrast_dir, "GSEA/Log2FC/MSIGDB/")
    dir.create(msigdb_dir, showWarnings = FALSE, recursive = TRUE)
    
    ib$visualizeGSEA(geneList = genelist, 
             go_class = "", 
             n_terms = 20, 
             outdir = msigdb_dir, 
             comp=comp,
             collection = "MSIGDB",
             org = "hsa",
             orgdb = org.Hs.eg.db,
             msigdb_file = "data/Human_MSigdb_March_01_2021_Entrezgene.gmt")
    
    
    # sign(Log2Fc)*-log10(Pvalue)
    #############################

    genelist <- sort(setNames(sign(res$log2FoldChange)*-log10(res$pvalue), rownames(res)), decreasing = TRUE)
    # genelist names should be entrezid
    names(genelist) <- mapIds(org.Hs.eg.db, 
                                keys = names(genelist), 
                                column = "ENTREZID", 
                                keytype = "ENSEMBL")
    genelist <- genelist[!is.na(names(genelist))]
    
    for(col in c("MF", "BP", "CC")){
      # GO BP Enrichment
      gsea_dir_tmp <- paste0(contrast_dir, "GSEA/signed_Pval/GO/", col, "/")
      dir.create(gsea_dir_tmp, showWarnings = FALSE, recursive = TRUE)
      
      ib$visualizeGSEA(geneList = genelist, 
               go_class = col, 
               n_terms = 20, 
               outdir = gsea_dir_tmp, 
               comp=comp,
               collection = "GO",
               org = "hsa",
               orgdb = org.Hs.eg.db,
               msigdb_file = "data/Human_MSigdb_March_01_2021_Entrezgene.gmt")
    }
    
    # KEGG Enrichment
    kegg_dir <- paste0(contrast_dir, "GSEA/signed_Pval/KEGG/")
    dir.create(kegg_dir, showWarnings = FALSE, recursive = TRUE)
    
    ib$visualizeGSEA(geneList = genelist, 
             go_class = "", 
             n_terms = 20, 
             outdir = kegg_dir, 
             comp=comp,
             collection = "KEGG",
             org = "hsa",
             orgdb = org.Hs.eg.db,
             msigdb_file = "data/Human_MSigdb_March_01_2021_Entrezgene.gmt")
    
    # MSIGDB Enrichment
    msigdb_dir <- paste0(contrast_dir, "GSEA/signed_Pval/MSIGDB/")
    dir.create(msigdb_dir, showWarnings = FALSE, recursive = TRUE)
    
    ib$visualizeGSEA(geneList = genelist, 
             go_class = "", 
             n_terms = 20, 
             outdir = msigdb_dir, 
             comp=comp,
             collection = "MSIGDB",
             org = "hsa",
             orgdb = org.Hs.eg.db,
             msigdb_file = "data/Human_MSigdb_March_01_2021_Entrezgene.gmt")
    
  }
}
