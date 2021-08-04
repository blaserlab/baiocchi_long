
top200 <- top200_by_specimen %>%
  filter(marker_test_q_value < 0.05) %>%
  group_by(gene_short_name) %>%
  summarise %>%
  pull

top200_by_specimen_matrix <- as.matrix(aggregate_gene_expression(cds_human_pass_counts[rowData(cds_human_pass_counts)$gene_short_name %in% top200,],
                          cell_group_df = data.frame(cell_id = rownames(colData(cds_human_pass_counts)),
                                                     specimen = colData(cds_human_pass_counts)$specimen)))
identical(rownames(top200_by_specimen_matrix), rownames(rowData(cds_human_pass_counts[rowData(cds_human_pass_counts)$gene_short_name %in% top200,])))#true
rownames(top200_by_specimen_matrix) <- rowData(cds_human_pass_counts[rowData(cds_human_pass_counts)$gene_short_name %in% top200,])$gene_short_name

top200_col_fun <- colorRamp2(breaks = c(min(scale(t(top200_by_specimen_matrix))),0,max(scale(t(top200_by_specimen_matrix)))), colors = heatmap_3_colors)

top200_highlights <- c("SRM", "CCND1", "PIK3CD", "MZB1", "TCL1A", "SUMO2", "GRB2")

top200_anno_df <-
  map(
    .x = top200_highlights,
    .f = function(x) {
      index <- which(colnames(t(top200_by_specimen_matrix)) == x)
      return(index)
    }
  ) %>% set_names(top200_highlights) %>%
  bind_cols() %>%
  pivot_longer(everything()) %>%
  as.data.frame()

top200_gene_anno <- HeatmapAnnotation(
  foo = anno_mark(
    at = top200_anno_df$value,
    labels = top200_anno_df$name,
    labels_gp = gpar(fontsize = 8),
    padding = 1.5,
    labels_rot = 45
  ),
  which = "column"
)

Heatmap(scale(t(top200_by_specimen_matrix)), name = "Top\nSpecific\nGenes",
        col = top200_col_fun,
        top_annotation = top200_gene_anno,
        column_dend_side = "bottom",
        show_column_names = F
        )

