use(ggplot2, utils, stats)

# Plotting Function GSEA
#' @export
fgsea_bars <- function(x,
                        select=6,
                        anot = comp) {
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