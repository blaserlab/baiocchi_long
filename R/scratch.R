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

