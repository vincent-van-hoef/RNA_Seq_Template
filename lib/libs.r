# Plotting Function GSEA
#' @export
fgsea_bars <- function(x,
                        select=6,
                        anot = comp) {

# Import
box::use(ggplot2, utils, stats)

  tmp <- x[order(x$NES), ]
  dtfr <- rbind(head(tmp, select), tail(tmp, select))
  dtfr$pathway <- gsub("%.*", "", dtfr$Description)
  dtfr$col <- ifelse(dtfr$NES < 0, "blue", "red")
  dtfr$col[dtfr$padj < 0.05 & dtfr$NES < 0] <- "blue4"
  dtfr$col[dtfr$padj < 0.05 & dtfr$NES > 0] <- "red4"
  dtfr$hjust <- ifelse(dtfr$NES > 0, 1, 0)
  ggplot(dtfr,
            aes(reorder(pathway, NES),
                NES,
                label = pathway,
                fill = col,
                hjust = hjust)) +
    geom_bar(stat = "identity") +
    scale_fill_identity() +
    geom_text(size = 3, aes(y = 0)) +
    coord_flip() +
    theme_classic() +
    theme(legend.position = "none", 
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),
          ends = "both"))) +
    ylab(paste(anot[3], "NES", anot[2], sep = "                                                     "))
}

#' @export
convertID <- function(df, type, keys, db) {

  # Import
  box::use(AnnotationDbi[mapIds]) 

  geneSymbols 	<- mapIds(db,
                            keys = rownames(df),
                            column = type,
                            keytype = keys,
                            multiVals = "first")
  inds 		<- which(!is.na(geneSymbols))
  found_genes 	<- geneSymbols[inds]
  df2 		<- df[names(found_genes), ]
  rownames(df2) 	<- found_genes
  return(df2)
}

# Create heatmap of vst or rld normalized dds object
#' @export
plotCustomHeatmap <- function(obj = rld, # normalized deseq2 object
                              plotGenes =  list, # single list of gene vectors,
                              anotationColumn = "Status",
                              anotationColor = list(Status = c("Control" = "green", "Fibro" = "purple")),
                              convertToSymbol = TRUE, # Convert to symbols first?
                              convertFromID = "ENSEMBL", # Convert from this ID
                              db = org.Hs.eg.db, # Change according to data, only for conversion
                              groupColumns = FALSE
                              ){

# Import
box::use(ComplexHeatmap[Heatmap, HeatmapAnnotation],
          SummarizedExperiment[assay, colData],
          grid[gpar],
          ./libs[convertID])

# Convert ID if necessary
if(convertToSymbol == TRUE) {
  convertedObj <- convertID(obj,
                            type = "SYMBOL",
                            keys = convertFromID,
                            db = db)
} else {
  convertedObj <- obj
  }

set <- names(plotGenes)
mat <- assay(convertedObj)[rownames(assay(convertedObj)) %in% plotGenes[[set]], ]
scaledMat <- t(scale(t(mat)))
ha = HeatmapAnnotation(Status = colData(obj)[, anotationColumn], col = anotationColor)
if (groupColumns == FALSE) {
  heatmapPlot <- Heatmap(scaledMat,
                       name = "Z-score", 
                       top_annotation = ha,
                       column_title = set,
                       row_names_gp = gpar(fontsize = 6),
                       column_names_gp = gpar(fontsize = 6),
                       rect_gp = gpar(col = "white", lwd = 1))
} else {
  heatmapPlot <- Heatmap(scaledMat,
                         name = "Z-score",
                         top_annotation = ha,
                         column_title = set,
                         row_names_gp = gpar(fontsize = 6),
                         column_names_gp = gpar(fontsize = 6),
                         rect_gp = gpar(col = "white", lwd = 1),
                         column_split = colData(obj)[, anotationColumn])
}
print(heatmapPlot)
}

# GO Enrichment
#' @export
gsea_viz <- function(geneList = genelist,
                    geneListID = "ENTREZID", # Should be ENTREZID
                     go_class = "BP",
                     n_terms = 20,
                     outdir = gobp_dir,
                     comp=comp,
                     collection = "GO", # One of c("GO", "KEGG", "MSIGDB")
                     org = "hsa",
                     msigdb_file = "data/Human_MSigdb_March_01_2021_Entrezgene.gmt"
                     ) {

# Import
box::use(enrichplot,
      clusterProfiler,
      org.Hs.eg.db,
      org.Mm.eg.db,
      org.Dr.eg.db,
      utils,
      grDevices,
      ./libs[fgsea_bars])

if (org == "hsa") {
    orgdb     <- org.Hs.eg.db
} else if (org == "mmu") {
   orgdb     <- org.Mm.eg.db
} else if (org == "dre") {
   orgdb     <- org.Dr.eg.db
}

 if (collection == "GO") {
  res_gsea <- gseGO(geneList     = geneList,
                     OrgDb        = orgdb,
                     keyType      = geneListID,
                     ont          = go_class,
                     minGSSize    = 10,
                     maxGSSize    = 500,
                     pvalueCutoff = 1,
                     verbose      = FALSE)
  } else if (collection == "KEGG") {
    res_gsea <- gseKEGG(geneList   = geneList,
                        organism     = org,
                        minGSSize    = 10,
                        maxGSSize    = 500,
                        pvalueCutoff = 1,
                        verbose      = FALSE)
  } else if (collection == "MSIGDB") {
    gmt <- clusterProfiler::read.gmt(gmtfile = msigdb_file)
    gmt$term <- gsub("%.*", "", gmt$term)
    res_gsea <- GSEA(geneList = geneList,
                     TERM2GENE = gmt,
                     pvalueCutoff = 1,
                     verbose = FALSE)
  }

res_gsea <- setReadable(res_gsea, orgdb, geneListID)
write.csv(res_gsea, paste0(paste0(outdir, paste(comp, collapse="_"), "_", collection, "_enrichment.csv")))
  
# Dotplot
pdf(paste0(outdir, paste(comp, collapse="_"), "_", collection, "_", go_class, "_dotplot.pdf"))
print(dotplot(res_gsea, showCategory=n_terms, font.size=10, color = "pvalue", x = "NES", label_format =))
dev.off()
  
# barplot
pdf(paste0(outdir, paste(comp, collapse="_"), "_", collection, "_", go_class, "_barplot.pdf"))
print(gsea_bars(x = res_gsea, select = 6, anot = comp))
dev.off()
  
# Enrichment map
emap <- enrichplot::pairwise_termsim(res_gsea, method = "JC", semData = NULL, showCategory = n_terms)
pdf(paste0(outdir, paste(comp, collapse="_"), "_", collection, "_", go_class, "_emap.pdf"))
print(emapplot(emap, showCategory=n_terms, layout = "kk", color = "pvalue"))
dev.off()
  
# GSEA plots
dir_gsea <- paste0(outdir, "GSEA_PLOTS/")
dir.create(dir_gsea, showWarnings = FALSE, recursive = TRUE)
for(i in 1:n_terms){
    pdf(paste0(dir_gsea, paste(comp, collapse="_"), "_", collection, "_", go_class, "_gseaplot_", gsub(" |/", "_", substr(res_gsea$Description[i], 0, 15)),  ".pdf"))
    print(gseaplot2(res_gsea, geneSetID = i, title = res_gsea$Description[i]))
    dev.off()
}
}

#' @export
plotVolcano <- function(res,
                        sig,
                        fc) {

# Import
box::use(ggplot2)

df <- as.data.frame(res)
df$id <- rownames(df)

# reduce too high significances to 10 on the Y-axis
df$pvalue <- ifelse(-log10(df$pvalue) > 10, 0.0000000001, df$pvalue)

#Number of significant genes
n_up <- nrow(subset(df, df$log2FoldChange > fc & df$pvalue < sig))
n_down <- nrow(subset(df, df$log2FoldChange < -fc & df$pvalue < sig))

df$col <- ifelse(abs(df$log2FoldChange) > fc & df$padj < sig,
                "Significant",
                "Not Significant")
df$col[is.na(df$col)] <- "Not Significant"

p1 <- ggplot(df, aes(log2FoldChange, -log10(pvalue))) +
    geom_point(size = 0.5, aes(col = col)) +
    scale_color_manual(values = c("black", "red")) +
    geom_hline(yintercept = -log10(sig), col = "grey") +
    geom_vline(xintercept = c(-fc,fc), col = "grey") +
    annotate("text", 
            x = max(df$log2FoldChange),
            y = -log10(sig / 2),
            label = n_up) +
    annotate("text", 
            x = min(df$log2FoldChange),
            y = -log10(sig / 2),
            label = n_down) +
    theme_classic() +
    theme(legend.position = "none") +
    ggtitle(gsub(".*:\\s[A-Za-z0-9]*\\s",
            "",
            res@elementMetadata@listData$description[2]))

p1
}