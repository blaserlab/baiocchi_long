# Figure 4A density umap faceted by treatment time -------------------------------
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

Four_A_density_umap_faceted <- bb_var_umap(cds_human_pass_sf,
                                    "density",
                                    facet_by = "treatment_time",
                                    sample_equally = F) +
  labs(color = "Cell\nDensity")


plot(Four_A_density_umap_faceted)

cowplot::save_plot(plot = Four_A_density_umap_faceted,
                   filename = fs::path(figs_out, "Four_A_density_umap_faceted", ext = "tiff"),
                   base_width = 10,
                   base_height = 4)

#supplemental fig 4A umap faceted by sample------------------------------------------
colData(cds_human_pass_sf)$treatment_time <-
  recode(
    colData(cds_human_pass_sf)$orig_id,
    "P5" = "Untreated Parental",
    "P61" = "Resistant PRMT5i",
    "P10" = "Short-term PRMT5i",
    "P53" = "Resistant PRMT5i",
    "P6_8_SP_VC_SURV" = "Short-Term Untreated",
    "P6_21_SP_VC" = "Long-term Untreated"
  )

Sup_Four_A_sample_umap_faceted <- bb_var_umap(cds_human_pass_sf, "sample",
             facet_by = c("treatment", "orig_id"), rows = vars(treatment), cols = vars(orig_id))



plot(Sup_Four_A_sample_umap_faceted)

cowplot::save_plot(plot = Sup_Four_A_sample_umap_faceted,
                   filename = fs::path(figs_out, "Sup_Four_A_sample_umap_faceted", ext = "tiff"),
                   base_width = 10,
                   base_height = 4)


# Figure 4B umap showing leiden clusters colored by enrichment ---------------------------
Four_B_leiden_enrichment_umap <- bb_var_umap(cds_human_pass_sf,
            "leiden",
            alt_label_col = "leiden",
            overwrite_labels = T, text_geom = "label") +
  labs(fill = "Log<sub>2</sub> Fold<br>Enrichment") +
  theme(legend.title = ggtext::element_markdown())

plot(Four_B_leiden_enrichment_umap)

cowplot::save_plot(plot = Four_B_leiden_enrichment_umap,
                   filename = fs::path(figs_out, "Four_B_leiden_umap", ext = "tiff"),
                   base_width = 3.5,
                   base_height = 3.25)

#Figure 4C Leiden PRMT5i Resistant vs. Short Term Enrichment--------------------------------------------

colData(cds_human_pass_sf)$treatment_2 <- recode(colData(cds_human_pass_sf)$sample,
                                                 "P1" = "Resistant PRMT5i",
                                                 "P2" = "Short-term PRMT5i",
                                                 "P3" = "Resistant PRMT5i",
                                                 "U1" = "Vehicle",
                                                 "U2" = "Vehicle",
                                                 "U3" = "Vehicle")

Four_C_leiden_PRMT5i_Resistant_vs_short_term_enrichment <-bb_cluster_representation2(
  obj = filter_cds(cds_human_pass_sf,
                   cells = bb_cellmeta(cds_human_pass_sf) |>
                     filter(treatment_2 %in% c("Resistant PRMT5i", "Short-term PRMT5i"))),
  sample_var = "sample",
  cluster_var = "leiden",
  comparison_var = "treatment_2",
  comparison_levels = c("Short-term PRMT5i", "Resistant PRMT5i"),
  sig_val = "FDR",
  return_val = "plot" # use "data" to get the table instead of the plot
)

plot(Four_C_leiden_PRMT5i_Resistant_vs_short_term_enrichment)

cowplot::save_plot(plot = Four_C_leiden_PRMT5i_Resistant_vs_short_term_enrichment,
                   filename = fs::path(figs_out, "Four_C_leiden_PRMT5i_Resistant_vs_short_term_enrichment", ext = "tiff"),
                   base_width = 6,
                   base_height = 4)

#Supplemental Figure 4B Leiden PRMT5i Resistant vs. Vehicle Enrichment--------------------------------------------
Sup_Four_B_leiden_PRMT5i_Resistant_vs_vehicle_enrichment <-bb_cluster_representation2(
  obj = filter_cds(cds_human_pass_sf,
                   cells = bb_cellmeta(cds_human_pass_sf) |>
                     filter(treatment_2 %in% c("Resistant PRMT5i", "Vehicle"))),
  sample_var = "sample",
  cluster_var = "leiden",
  comparison_var = "treatment_2",
  comparison_levels = c("Vehicle", "Resistant PRMT5i"),
  sig_val = "FDR",
  return_val = "plot" # use "data" to get the table instead of the plot
)


plot(Sup_Four_B_leiden_PRMT5i_Resistant_vs_vehicle_enrichment)

cowplot::save_plot(plot = Sup_Four_B_leiden_PRMT5i_Resistant_vs_vehicle_enrichment,
                   filename = fs::path(figs_out, "Sup_Four_B_leiden_PRMT5i_Resistant_vs_vehicle_enrichment", ext = "tiff"),
                   base_width = 6,
                   base_height = 4)

# Figure 4D gene set heat map -----------------------------------------------

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

heatmap_colfun <-
  circlize::colorRamp2(breaks =  c(min(agg_mat), 0, max(agg_mat)), colors = heatmap_3_colors)

Four_D_pathway_heatmap <- grid.grabExpr(draw(
  Heatmap(
    t(agg_mat),
    name = "Pathway\nExpression",
    column_names_rot = 30,
    column_dend_height = unit(2, "mm"),
    row_dend_width = unit(2, "mm"),
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    # right_annotation = heatmap_rowanno,
    row_title = "Pathway",
    column_title = "Leiden Cluster",
    column_title_side = "bottom",
    col = heatmap_colfun,
    heatmap_legend_param = list(
      direction = "horizontal",
      row_title_position = "left", column_title_position = "bottom")
  ), heatmap_legend_side = "bottom"
), wrap = T)

plot_grid(Four_D_pathway_heatmap)

cowplot::save_plot(plot = Four_D_pathway_heatmap,
                   filename = fs::path(figs_out, "Four_D_pathway_heatmap", ext = "tiff"),
                   base_width = 4,
                   base_height = 2)


#Fig 4E _mtor_upregulation_clusters

Four_E_mtor_umap <- bb_gene_umap(cds_human_pass_sf, gene_or_genes = "MTOR")


plot_grid(Four_E_mtor_umap)

cowplot::save_plot(plot = Four_E_mtor_umap,
                   filename = fs::path(figs_out, "Four_E_mtor_umap", ext = "tiff"),
                   base_width = 6,
                   base_height = 4)

