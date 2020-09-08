# Functions for plotting files:
# 1. PCA_plot
# 2. volcano
# 3. gene_counts


####################################################
# 1. PCA_plot (plots PCA results, colors points according to given category)
PCA_plot <- function(pcaData, category) {
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color=pcaData[,category])) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    theme(legend.title = element_blank())
}


####################################################
# 2. volcano (creates volcano plot that labels top 10 genes)
volcano <- function(res_table_tb, padj_cutoff, lfc_cutoff) {
  # Add a column indicating whether each gene is significant; order the rows by padj; create new "genelabels" column.
  res_table_tb <- res_table_tb %>%
    mutate(threshold = padj < padj_cutoff & abs(log2FoldChange) >= lfc_cutoff) %>%
    arrange(padj) %>%
    mutate(genelabels = "")
  # Fill in the genelabels column with gene names for the top genes
  res_table_tb$genelabels[1:10] <- res_table_tb$gene_name[1:10]
  # Plot
  ggplot(res_table_tb, aes(x=log2FoldChange, y=-log10(padj))) +
    geom_point(aes(colour=threshold)) +
    scale_color_manual(values=c("black", "red")) +
    geom_text_repel(aes(label=genelabels)) +
    xlab("log2 fold change") +
    ylab("-log10 adjusted p-value") +
    theme(legend.position = "none", plot.title=element_text(size=rel(1.5), hjust=0.5), axis.title=element_text(size=rel(1.25)))
}


####################################################
# 3. gene_counts (plots small RNA levels for a specific gene)
gene_counts <- function(dds, gene, sample_order) {
  d <- plotCounts(dds, gene = gene, intgroup = c("Genotype", "Replicate"), returnData = TRUE)
  d$Replicate <- as.character(d$Replicate)
  d$Genotype <- factor(d$Genotype, levels = sample_order)
  gene_name <- gene_name <- gene_names[which(gene_names$gene == gene), 2]
  ggplot(d, aes(fill = Replicate, x = Genotype, y = count)) +
    geom_bar(position = "dodge", stat = "identity") +
    theme_bw() +
    ggtitle(gene_name) +
    ylab("normalized counts") +
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
}
