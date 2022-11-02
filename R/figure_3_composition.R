# figure 3 ----------------------------------------
#mackenzie_comment
fig_3_top <- plot_grid(NULL)

fig_3_mid <- plot_grid(
  ggrastr::rasterise(density_umap_faceted, dpi = 300),
  ncol = 1,
  labels = c("D)")
)

fig_3_bottom_right <- plot_grid(
  leiden_enrichment_barplot +
    coord_flip() +
    theme(axis.title.x = ggtext::element_markdown()),
  pathway_heatmap,
  ncol = 1,
  rel_heights = c(1,1),
  labels = c("F)", "G)"), hjust = 1.5

)

fig_3_bottom <- plot_grid(
  ggrastr::rasterise(leiden_enrichment_umap, dpi = 300),
  fig_3_bottom_right,
  ncol = 2,
  rel_widths = c(1.2, 1),
  labels = c("E)", "")
)

fig_3 <- plot_grid(
  fig_3_top,
  fig_3_mid,
  fig_3_bottom,
  nrow = 3,
  rel_heights = c(1.2, 1, 1),
  align = "v",
  axis = "l"
)

cowplot::save_plot(plot = fig_3,
                   filename = fs::path(figs_out, "figure_3", ext = "pdf"),
                   base_width = 7.5,
                   base_height = 9.75)
