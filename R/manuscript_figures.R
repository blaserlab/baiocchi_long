# density umap faceted by treatment time -------------------------------
colData(cds_human_pass_sf)$treatment_time <-
  recode(
    colData(cds_human_pass_sf)$orig_id,
    "P5" = "Untreated",
    "P61" = "Long term treated",
    "P10" = "Short term treated",
    "P53" = "Long term treated",
    "P6_8_SP_VC_SURV" = "Untreated",
    "P6_21_SP_VC" = "Untreated"
  )

colData(cds_human_pass_sf)$treatment_time <-
  factor(colData(cds_human_pass_sf)$treatment_time,
         levels = c("Untreated", "Short term treated", "Long term treated"))

density_umap_faceted <- bb_var_umap(cds_human_pass_sf,
                                    "density",
                                    facet_by = "treatment_time",
                                    sample_equally = F) +
  labs(color = "Cell\nDensity")

# umap showing leiden clusters colored by enrichment ---------------------------
# calculate the enrichment
leiden_long_short_enrichment_tbl <- bb_cluster_representation(
  cds = cds_human_pass_sf[,colData(cds_human_pass_sf)$treatment_time %in% c("Long term treated", "Short term treated")],
  cluster_var = "leiden",
  class_var = "treatment_time",
  experimental_class = "Long term treated",
  control_class = "Short term treated",
  return_value = "table"
)

leiden_long_short_enrichment <- leiden_long_short_enrichment_tbl |>
  select(leiden, log2fold_change_over_control) |>
  deframe()

# pass that into the cds as a new metadata column
colData(cds_human_pass_sf)$leiden_long_short_enrichment <- recode(colData(cds_human_pass_sf)$leiden, !!!leiden_long_short_enrichment)

leiden_enrichment_umap <- bb_var_umap(cds_human_pass_sf,
            "leiden_long_short_enrichment",
            alt_label_col = "leiden",
            overwrite_labels = T, text_geom = "label") +
  theme(legend.position = "right") +
  labs(fill = "Log<sub>2</sub> Fold<br>Enrichment") +
  theme(legend.title = ggtext::element_markdown())

leiden_enrichment_barplot <- ggplot(leiden_long_short_enrichment_tbl,
       mapping = aes(x = reorder(leiden, desc(log2fold_change_over_control)),
                     y = log2fold_change_over_control,
                     fill = log2fold_change_over_control)) +
  geom_col(color = "black") +
  scale_fill_viridis_c() +
  theme(legend.position = "none") +
  labs(x = "Leiden Cluster", y = "Log<sub>2</sub> Fold Enrichment") +
  theme(axis.title.y = ggtext::element_markdown()) +
  geom_text(mapping = aes(y = texty, label = p.signif), size = 3, show.legend = F, vjust = -0.5) +
  expand_limits(y = 6)


# gene set umaps ----------------------------------------------------
bb_gene_umap(cds_human_pass_sf, bb_rowmeta(cds_human_pass_sf) |> select(feature_id, IGF1_signaling) |> filter(IGF1_signaling), max_expr_val = 3.5)


#figure_3G
# # add gene sets provided by mackenzie
# mac_genesets_to_add <- map2_dfr(
#   .x = list.files(
#     "~/network/X/Labs/Blaser/share/collaborators/baiocchi_long_manuscript/genesets/xls3",
#     full.names = TRUE
#   ),
#   .y = c(
#     "mtor",
#     "igf",
#     "regulation_of_eif4_and_p70S6k",
#     "estrogen_receptor"
#   ),
#   .f = \(x, y, cds = cds_human_pass_sf) {
#     readxl::read_excel(x) |>
#       left_join(bb_rowmeta(cds),
#                 by = c("Symbol" = "gene_short_name")) |>
#       filter(!is.na(feature_id)) |>
#       select(feature_id) |>
#       mutate(mac_geneset = TRUE) |>
#       right_join(bb_rowmeta(cds)) |>
#       select(feature_id, mac_geneset) |>
#       pivot_longer(mac_geneset) |>
#       mutate(value = ifelse(is.na(value), FALSE, TRUE)) |>
#       mutate(name = y)
#
#
#   }
# ) |> pivot_wider(names_from = "name", values_from = "value")
#
# cds_human_pass_sf <- bb_tbl_to_rowdata(obj = cds_human_pass_sf, min_tbl = mac_genesets_to_add)
#
# dir.create("data")
#
# save(cds_human_pass_sf, file = "data/cds_human_pass_sf.rda", compress = "bzip2")
#
# bb_rowmeta(cds_human_pass_sf)
#
# # produce the plots
# walk2(
#   .x = c(
#     "mtor",
#     "igf",
#     "regulation_of_eif4_and_p70S6k",
#     "estrogen_receptor"
#
#   ),
#   .y = c(Inf, 4, 6, 10, Inf, 3.5, Inf, 2.5, 4, 3),
#   .f = \(x, y, dat = cds_human_pass_sf, figdir = figs_out) {
#     p <-
#       bb_gene_umap(
#         dat,
#         bb_rowmeta(dat) |> filter(!!sym(x) &
#                                     !is.na(module))
#         |> select(feature_id, !!sym(x)),
#         max_expr_val = y
#       ) +
#       theme(strip.text = element_blank()) +
#       labs(title = x)
#     p <- ggrastr::rasterise(p, dpi = 300)
#     save_plot(
#       filename = str_glue("{figdir}/mac_geneset_figs/{x}.pdf"),
#       plot = p,
#       base_width = 5.4,
#       base_height = 4.2
#     )
#   }
# )
#
# # print out the gene sets used
# bb_rowmeta(cds_human_pass_sf) |>
#   select(
#     gene_short_name,
#     mtor,
#     igf,
#     regulation_of_eif4_and_p70S6k,
#     estrogen_receptor
#   ) |>
#   pivot_longer(-gene_short_name, names_to = "pathway") |>
#   filter(value) |>
#   select(gene_short_name, pathway) |>
#   arrange(pathway) |>
#   write_csv(str_glue("{figs_out}/mac_geneset_figs/pathways.csv"))


agg_mat_list <- map(
  .x = c(
    "mtor",
    "IGF1_signaling",
    "eif4_p70S6k",
    "ER_signaling"

  ),
  .f = \(x, dat = cds_human_pass_sf) {
    mat <- aggregate_gene_expression(
      dat,
      gene_group_df = bb_rowmeta(dat) |> select(feature_id, !!sym(x)),
      cell_group_df = bb_cellmeta(dat) |> select(cell_id, leiden),
      scale_agg_values = F
    ) |>
      t()
    mat <- matrix(mat[, "TRUE"])
    colnames(mat) <- x
    return(mat)

  }
)

agg_mat <- bind_cols(agg_mat_list) |>
  scale()
rownames(agg_mat) <- paste0("leiden ", 1:13)

heatmap_rowanno_df <- bb_cellmeta(cds_human_pass_sf) |>
  group_by(leiden, partition) |>
  summarise() |>
  mutate(leiden = paste0("leiden ", leiden)) |>
  column_to_rownames("leiden")

heatmap_colfun <-
  circlize::colorRamp2(breaks =  c(min(agg_mat), 0, max(agg_mat)), colors = heatmap_3_colors)

pathway_heatmap <- grid.grabExpr(draw(
  Heatmap(
    agg_mat,
    name = "Pathway\nExpression",
    column_names_rot = 30,
    # right_annotation = heatmap_rowanno,
    col = heatmap_colfun,
    column_title = "Pathway",
    row_title = "Cell Cluster"
  )
), wrap = T)

plot_grid(pathway_heatmap)

save_plot(
  pathway_heatmap,
  filename = str_glue("{figs_out}/mac_geneset_figs/pathway_heatmap.pdf"),
  base_height = 5,
  base_width = 5
)
