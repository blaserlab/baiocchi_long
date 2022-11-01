#figure_3D
analysis_configs
colData(cds_human_pass_sf)$treatment_time <-
  recode(
    colData(cds_human_pass_sf)$orig_id,
    "P5" = "Untreated",
    "P61" = "long_term_treated",
    "P10" = "short_term_treated",
    "P53" = "long_term_treated",
    "P6_8_SP_VC_SURV" = "Untreated",
    "P6_21_SP_VC" = "Untreated"
  )
colData(cds_human_pass_sf)$treatment_time <- factor(colData(cds_human_pass_sf)$treatment_time, levels = c("Untreated", "short_term_treated", "long_term_treated"))
bb_var_umap(cds_human_pass_sf, "density", facet_by = "treatment_time", sample_equally = F)

#figure_3E
bb_var_umap(cds_human_pass_sf, "leiden")

#figure_3F
analysis_configs
colData(cds_human_pass_sf)$treatment_time <-
  recode(
    colData(cds_human_pass_sf)$orig_id,
    "P5" = "Untreated",
    "P61" = "long_term_treated",
    "P10" = "short_term_treated",
    "P53" = "long_term_treated",
    "P6_8_SP_VC_SURV" = "Untreated",
    "P6_21_SP_VC" = "Untreated"
  )
bb_cluster_representation(
  cds = cds_human_pass_sf[,colData(cds_human_pass_sf)$treatment_time %in% c("long_term_treated", "short_term_treated")],
  cluster_var = "leiden",
  class_var = "treatment_time",
  experimental_class = "long_term_treated",
  control_class = "short_term_treated",
  return_value = "table"
)

bb_cluster_representation(
  cds = cds_human_pass_sf[,colData(cds_human_pass_sf)$treatment_time %in% c("long_term_treated", "short_term_treated")],
  cluster_var = "leiden",
  class_var = "treatment_time",
  experimental_class = "long_term_treated",
  control_class = "short_term_treated",
  return_value = "plot"
)

#figure_3G
# add gene sets provided by mackenzie
mac_genesets_to_add <- map2_dfr(
  .x = list.files(
    "~/network/X/Labs/Blaser/share/collaborators/baiocchi_long_manuscript/genesets/xls3",
    full.names = TRUE
  ),
  .y = c(
    "mtor",
    "igf",
    "regulation_of_eif4_and_p70S6k",
    "estrogen_receptor"
  ),
  .f = \(x, y, cds = cds_human_pass_sf) {
    readxl::read_excel(x) |>
      left_join(bb_rowmeta(cds),
                by = c("Symbol" = "gene_short_name")) |>
      filter(!is.na(feature_id)) |>
      select(feature_id) |>
      mutate(mac_geneset = TRUE) |>
      right_join(bb_rowmeta(cds)) |>
      select(feature_id, mac_geneset) |>
      pivot_longer(mac_geneset) |>
      mutate(value = ifelse(is.na(value), FALSE, TRUE)) |>
      mutate(name = y)


  }
) |> pivot_wider(names_from = "name", values_from = "value")

cds_human_pass_sf <- bb_tbl_to_rowdata(obj = cds_human_pass_sf, min_tbl = mac_genesets_to_add)

dir.create("data")

save(cds_human_pass_sf, file = "data/cds_human_pass_sf.rda", compress = "bzip2")


# produce the plots
walk2(
  .x = c(
    "mtor",
    "igf",
    "regulation_of_eif4_and_p70S6k",
    "estrogen_receptor"

  ),
  .y = c(Inf, 4, 6, 10, Inf, 3.5, Inf, 2.5, 4, 3),
  .f = \(x, y, dat = cds_human_pass_sf, figdir = figs_out) {
    p <-
      bb_gene_umap(
        dat,
        bb_rowmeta(dat) |> filter(!!sym(x) &
                                    !is.na(module))
        |> select(feature_id, !!sym(x)),
        max_expr_val = y
      ) +
      theme(strip.text = element_blank()) +
      labs(title = x)
    p <- ggrastr::rasterise(p, dpi = 300)
    save_plot(
      filename = str_glue("{figdir}/mac_geneset_figs/{x}.pdf"),
      plot = p,
      base_width = 5.4,
      base_height = 4.2
    )
  }
)

# print out the gene sets used
bb_rowmeta(cds_human_pass_sf) |>
  select(
    gene_short_name,
    mtor,
    igf,
    regulation_of_eif4_and_p70S6k,
    estrogen_receptor
  ) |>
  pivot_longer(-gene_short_name, names_to = "pathway") |>
  filter(value) |>
  select(gene_short_name, pathway) |>
  arrange(pathway) |>
  write_csv(str_glue("{figs_out}/mac_geneset_figs/pathways.csv"))


agg_mat_list <- map(
  .x = c(
    "mtor",
    "igf",
    "regulation_of_eif4_and_p70S6k",
    "estrogen_receptor"

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
    right_annotation = heatmap_rowanno,
    col = heatmap_colfun,
    column_title = "Pathway",
    row_title = "Cell Cluster"
  )
), wrap = T)


save_plot(
  pathway_heatmap,
  filename = str_glue("{figs_out}/mac_geneset_figs/pathway_heatmap.pdf"),
  base_height = 5,
  base_width = 5
)
