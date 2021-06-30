cds_human_top_markers %>%
  arrange(cluster_method, cell_group) %>%
  write_csv(str_glue("{tables_out}/cds_human_top_markers.csv"))
