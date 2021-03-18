import("ComplexHeatmap", "modules")
export("plotCustomHeatmap")

# Create heatmap of vst or rld normalized dds object
plotCustomHeatmap <- function(obj = rld, # normalized deseq2 object
                              plotGenes =  list, # single list of gene vectors,
                              anotationColumn = "Status",
                              anotationColor = list(Status = c("Control" = "green", "Fibro" = "purple")),
                              convertToSymbol = TRUE, # Convert to symbols first?
                              convertFromID = "ENSEMBL", # Convert from this ID
                              db = org.Hs.eg.db, # Change according to data, only for conversion
                              groupColumns = FALSE
                              ){

# Convert ID if necessary
if(convertToSymbol == TRUE) {
  lib     <- modules::use("convertID.R")
  convertedObj <- lib$convertID(obj,
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