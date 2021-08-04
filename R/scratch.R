bb_gene_umap(cds_human_pass_counts, c("LYN","PIK3CD","MRTO4", "MYC", "SRM", "PLK1", "HSPE1", "MAD2L1", "KPNA2"))
plot_genes_in_pseudotime(cds_human_pass_counts[rowData(cds_human_pass_counts)$gene_short_name == "PIK3CD",])

bb_var_umap(cds_human_pass_counts, "pseudotime")
bb_var_umap(cds_human_pass_counts, "specimen")
bb_var_umap(cds_human_pass_counts, "leiden")
cds_human_pass_counts_top_markers %>% filter(cell_group == "leiden 4") %>% View()



# if you want to plot aggregate expression of groups of genes
# create a data frame with columns for gene_id and gene_group
# then feed it into bb_gene_umap.  For example:

bb_gene_umap(
  cds = cds_human_pass_counts,
  gene_or_genes = data.frame(gene_id = rowData(cds_human_pass_counts)$id,
                             gene_grouping = rowData(cds_human_pass_counts)$module_labeled)
)

# you'll have to use some of your coding skills to translate from gene_short_name to gene_id
# and get your data table in this form.
# bb_gene_umap will produce a faceted plot with a facet for each grouping
# if you have one or more solid lists of gene groupings, like say from the literature
# or from some other experiment we can add them into the gene metadata table
# and make your code a little cleaner and error-proof
