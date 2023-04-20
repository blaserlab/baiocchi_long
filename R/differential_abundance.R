colData(cds_human_pass_sf)$treatment_1 <- recode(colData(cds_human_pass_sf)$sample,
                                                 "P1" = "long-term",
                                                 "P2" = "short-term/untreated",
                                                 "P3" = "long-term",
                                                 "U1" = "short-term/untreated",
                                                 "U2" = "short-term/untreated",
                                                 "U3" = "short-term/untreated")

colData(cds_human_pass_sf)$treatment_2 <- recode(colData(cds_human_pass_sf)$sample,
                                                 "P1" = "long-term",
                                                 "P2" = "short-term",
                                                 "P3" = "long-term",
                                                 "U1" = "untreated",
                                                 "U2" = "untreated",
                                                 "U3" = "untreated")

# Differential abundance  http://bioconductor.org/books/3.13/OSCA.multisample/differential-abundance.html

# long-term vs short-term AND untreated together
bb_cluster_representation2(
  obj = cds_human_pass_sf,
  sample_var = "sample",
  cluster_var = "leiden",
  comparison_var = "treatment_1",
  comparison_levels = c("short-term/untreated", "long-term"),
  sig_val = "FDR",
  return_val = "plot" # use "data" to get the table instead of the plot
)

# long-term vs untreated ONLY.
# Short-term treated sample filtered out
bb_cluster_representation2(
  obj = filter_cds(cds_human_pass_sf,
                   cells = bb_cellmeta(cds_human_pass_sf) |>
                     filter(treatment_2 %in% c("long-term", "untreated"))),
  sample_var = "sample",
  cluster_var = "leiden",
  comparison_var = "treatment_2",
  comparison_levels = c("untreated", "long-term"),
  sig_val = "FDR",
  return_val = "plot" # use "data" to get the table instead of the plot
)

# long-term vs short_term treated ONLY.
# untreated samples filtered out
bb_cluster_representation2(
  obj = filter_cds(cds_human_pass_sf,
                   cells = bb_cellmeta(cds_human_pass_sf) |>
                     filter(treatment_2 %in% c("long-term", "short-term"))),
  sample_var = "sample",
  cluster_var = "leiden",
  comparison_var = "treatment_2",
  comparison_levels = c("short-term", "long-term"),
  sig_val = "FDR",
  return_val = "plot" # use "data" to get the table instead of the plot
)

