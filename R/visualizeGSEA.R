import("enrichplot", "clusterProfiler", "org.Hs.eg.db", "org.Mm.eg.db", "org.Dr.eg.db")
export("gsea_viz")

# GO Enrichment
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

lib     <- modules::use("R/barplotGsea.R")

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
    res_gsea <- GSEA(geneList = genelist,
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
print(lib$barplotGsea$fgsea_bars(x = res_gsea, select = 6, anot = comp))
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