import("enrichplot")
export("gsea_viz")

# GO Enrichment
gsea_viz <- function(geneList = genelist, 
                     go_class = "BP", 
                     n_terms = 20, 
                     outdir = gobp_dir, 
                     comp=comp,
                     collection = "GO",
                     org = "hsa",
                     orgdb = org.Hs.eg.db,
                     msigdb_file = "data/Human_MSigdb_March_01_2021_Entrezgene.gmt"
                     ) {
  
  if(collection == "GO"){
  res_gsea <- gseGO(geneList     = geneList,
                     OrgDb        = orgdb,
                     keyType      = "ENTREZID", 
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

  res_gsea <- setReadable(res_gsea, orgdb, 'ENTREZID')
  write.csv(res_gsea, paste0(paste0(outdir, paste(comp, collapse="_"), "_", collection, "_enrichment.csv")))
  
  # Dotplot
  pdf(paste0(outdir, paste(comp, collapse="_"), "_", collection, "_", go_class, "_dotplot.pdf"))
  print(dotplot(res_gsea, showCategory=n_terms, font.size=10, color = "pvalue", x = "NES", label_format =))
  dev.off()
  
  # barplot
  pdf(paste0(outdir, paste(comp, collapse="_"), "_", collection, "_", go_class, "_barplot.pdf"))
  print(fgsea_bars(x = res_gsea, select = 6, anot = comp))
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

getMatrixWithSelectedIds <- function(df, type, keys, db){
  geneSymbols 	<- mapIds(db, keys=rownames(df), column=type, keytype=keys, multiVals="first")
  inds 		<- which(!is.na(geneSymbols))
  found_genes 	<- geneSymbols[inds]
  df2 		<- df[names(found_genes), ]
  rownames(df2) 	<- found_genes
  return(df2)
}

