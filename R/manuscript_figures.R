# density umap faceted by treatment time -------------------------------
colData(cds_human_pass_sf)$treatment_time <-
  recode(
    colData(cds_human_pass_sf)$orig_id,
    "P5" = "Vehicle",
    "P61" = "Resistant PRMT5i",
    "P10" = "Short-term PRMT5i",
    "P53" = "Resistant PRMT5i",
    "P6_8_SP_VC_SURV" = "Vehicle",
    "P6_21_SP_VC" = "Vehicle"
  )

colData(cds_human_pass_sf)$treatment_time <-
  factor(colData(cds_human_pass_sf)$treatment_time,
         levels = c("Vehicle", "Short-term PRMT5i", "Resistant PRMT5i"))

density_umap_faceted <- bb_var_umap(cds_human_pass_sf,
                                    "density",
                                    facet_by = "treatment_time",
                                    sample_equally = F) +
  labs(color = "Cell\nDensity")



plot(density_umap_faceted)

cowplot::save_plot(plot = density_umap_faceted,
                   filename = fs::path(figs_out, "density_umap_faceted", ext = "tiff"),
                   base_width = 10,
                   base_height = 4)

# umap showing leiden clusters colored by enrichment ---------------------------
# calculate the enrichment
leiden_long_short_enrichment_tbl <- bb_cluster_representation(
  cds = cds_human_pass_sf[,colData(cds_human_pass_sf)$treatment_time %in% c("Resistant PRMT5i", "Short-term PRMT5i")],
  cluster_var = "leiden",
  class_var = "treatment_time",
  experimental_class = "Resistant PRMT5i",
  control_class = "Short-term PRMT5i",
  return_value = "table"
)

leiden_long_short_enrichment <- leiden_long_short_enrichment_tbl |>
  select(leiden, log2fold_change_over_control) |>
  deframe()

# pass that into the cds as a new metadata column
colData(cds_human_pass_sf)$leiden_long_short_enrichment <- recode(colData(cds_human_pass_sf)$leiden, !!!leiden_long_short_enrichment)

leiden_enrichment_umap <- bb_var_umap(cds_human_pass_sf,
            "leiden",
            alt_label_col = "leiden",
            overwrite_labels = T, text_geom = "label") +
  labs(fill = "Log<sub>2</sub> Fold<br>Enrichment") +
  theme(legend.title = ggtext::element_markdown())

plot(leiden_enrichment_umap)

cowplot::save_plot(plot = leiden_enrichment_umap,
                   filename = fs::path(figs_out, "leiden_umap", ext = "tiff"),
                   base_width = 3.5,
                   base_height = 3.25)

leiden_enrichment_barplot <- ggplot(leiden_long_short_enrichment_tbl,
       mapping = aes(x = reorder(leiden, desc(log2fold_change_over_control)),
                     y = log2fold_change_over_control,
                     fill = log2fold_change_over_control)) +
  geom_col(color = "black") +
  scale_fill_viridis_c() +
  theme(legend.position = "none") +
  labs(x = "Leiden Cluster", y = "Log<sub>2</sub> Fold Resistant vs. Short-term PRMT5i") +
  theme(axis.title.y = ggtext::element_markdown()) +
  geom_text(mapping = aes(y = texty, label = p.signif), size = 10, show.legend = F, vjust = -0.5) +
  expand_limits(y = 6)

plot(leiden_enrichment_barplot)

cowplot::save_plot(plot = leiden_enrichment_barplot,
                   filename = fs::path(figs_out, "leiden_enrichment_barplot", ext = "tiff"),
                   base_width = 15,
                   base_height = 10)
# gene set heat map -----------------------------------------------

agg_mat_list <- map(
  .x = c(
    "mtor",
    "IGF1_signaling",
    "eif4_p70S6k",
    "pi3k_akt", "ERK_MAPK_signaling"

  ),
  .f = \(x, dat = cds_human_pass_sf) {
    mat <- bb_aggregate(
      dat,
      gene_group_df = bb_rowmeta(dat) |> select(feature_id, !!sym(x)),
      cell_group_df = bb_cellmeta(dat) |> select(cell_id, leiden),
      scale_agg_values = F
    ) |>
      t()
    true_mat <- matrix(mat[, "TRUE"])
    colnames(true_mat) <- x
    rownames(true_mat) <- rownames(mat)
    return(true_mat)

  }
)

agg_mat <- bind_cols(agg_mat_list) |>
  scale()
rownames(agg_mat) <- paste0(1:13)

# heatmap_rowanno_df <- bb_cellmeta(cds_human_pass_sf) |>
#   group_by(leiden, partition) |>
#   summarise() |>
#   mutate(leiden = paste0("leiden ", leiden)) |>
#   column_to_rownames("leiden")

heatmap_colfun <-
  circlize::colorRamp2(breaks =  c(min(agg_mat), 0, max(agg_mat)), colors = heatmap_3_colors)

pathway_heatmap <- grid.grabExpr(draw(
  Heatmap(
    t(agg_mat),
    name = "Pathway\nExpression",
    column_names_rot = 30,
    column_dend_height = unit(2, "mm"),
    row_dend_width = unit(2, "mm"),
    row_names_gp = gpar(fontsize = 9),
    column_names_gp = gpar(fontsize = 9),
    # right_annotation = heatmap_rowanno,
    col = heatmap_colfun,
    heatmap_legend_param = list(
      direction = "horizontal",
      title_position = "lefttop")
  ), heatmap_legend_side = "bottom"
), wrap = T)

plot_grid(pathway_heatmap)

cowplot::save_plot(plot = pathway_heatmap,
                   filename = fs::path(figs_out, "pathway_heatmap", ext = "tiff"),
                   base_width = 4,
                   base_height = 1.75)

#subfigure_mtor_upregulation_clusters

bb_gene_umap(cds_human_pass_sf, gene_or_genes = "MTOR")

bb_gene_dotplot(cds_human_pass_sf[,colData(cds_human_pass_sf)$leiden %in% c("3","13","2","6")], markers = c("MTOR"), group_cells_by = "leiden")
