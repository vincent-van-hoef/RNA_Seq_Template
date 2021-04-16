box::use(ggplot2)

#' @export
plotVolcano <- function(res,
                        sig,
                        fc) {
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