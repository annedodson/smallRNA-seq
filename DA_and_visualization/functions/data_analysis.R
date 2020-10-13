# Functions for data analysis:
# 1. new_res_table
# 2. sig_genes


####################################################
# 1. new_res_table (converts results table to a tibble and add a column with meaningul gene names)
new_res_table <- function(res_table) {
  res_table_tb <- res_table %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble()
  res_table_tb <- inner_join(gene_names, res_table_tb)
  return(res_table_tb)
}


####################################################
# 2. sig_genes (creates a list of significant genes using specified cutoffs)
sig_genes <- function(res_table_tb, padj_cutoff, lfc_cutoff) {
  sig_tb <- res_table_tb %>% 
    filter(padj < padj_cutoff & abs(log2FoldChange) > lfc_cutoff)
  return(sig_tb)
}
