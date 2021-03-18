import{"AnnotationDbi"}
export("convertID")

convertID <- function(df, type, keys, db) {

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