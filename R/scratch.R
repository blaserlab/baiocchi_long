bb_var_umap(cds_human, "classification")
bb_var_umap(cds_human, "specimen")
bb_var_umap(cds_human, "Size_Factor")
bb_var_umap(cds_human, "partition")
bb_gene_umap(cds_human, "CD19")
bb_gene_umap(cds_human, "CCND1")

cds_human_top_markers %>% filter(cluster_method == "partition") %>% View()
colData(cds_human) %>%
  as_tibble() %>%
  group_by(specimen, classification, passage) %>%
  summarise(n = n())

bb_var_umap(cds = cds_human, var = "specimen")
bb_var_umap(cds = cds_human, var = "date_collected")
bb_var_umap(cds = cds_human, var = "specimen", alt_dim_x = "prealignment_dim1", alt_dim_y = "prealignment_dim2")
bb_var_umap(cds = cds_human, var = "date_collected", alt_dim_x = "prealignment_dim1", alt_dim_y = "prealignment_dim2")

analysis_configs

bb_var_umap(cds_human, "classification", foreground_alpha = 0.4)
bb_var_umap(cds_human, "pseudotime", foreground_alpha = 0.4)
colData(cds_human)

plot_cells(cds_human, color_cells_by = "pseudotime")

bb_var_umap(cds_human, "clonotype", facet_by = "value")
bb_var_umap(cds_human, "specimen", facet_by = "value")
bb_var_umap(cds_human, "clonotype", facet_by = "specimen")
bb_var_umap(cds_human, "date_collected", facet_by = "specimen")
bb_var_umap(cds_human, "partition")


cds_human_top_markers %>%
  filter(cluster_method == "partition") %>%
  filter(cell_group %in% c("partition 1", "partition 2")) %>%
  View()

bb_gene_umap(cds_human, "HNRNPH1")

bb_gene_dotplot(cds_human[,colData(cds_human)$partition %in% c("1", "2")], markers = c("HNRNPH1","CCND1", "MDM2", "MDM4"), group_cells_by = "partition")
bb_var_umap(cds_human, "leiden", overwrite_labels = T)

cds_human_top_markers %>%
  filter(cluster_method == "leiden") %>%
  View()

bb_gene_umap(cds_human, gene_or_genes = c("CD79A", "CD79B", "PIK3CD", "PRDM2"))
plot_genes_in_pseudotime(cds_subset = cds_human[rowData(cds_human)$gene_short_name %in% c("PIK3CD"),
                                                colData(cds_human)$partition == "2"])
?bb_var_umap
