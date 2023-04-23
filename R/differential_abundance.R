colData(cds_human_pass_sf)$treatment_1 <- recode(colData(cds_human_pass_sf)$sample,
                                                 "P1" = "Resistant PRMT5i",
                                                 "P2" = "Short-term PRMT5i/Vehicle",
                                                 "P3" = "Resistant PRMT5i",
                                                 "U1" = "Short-term PRMT5i/Vehicle",
                                                 "U2" = "Short-term PRMT5i/Vehicle",
                                                 "U3" = "Short-term PRMT5i/Vehicle")

colData(cds_human_pass_sf)$treatment_2 <- recode(colData(cds_human_pass_sf)$sample,
                                                 "P1" = "Resistant PRMT5i",
                                                 "P2" = "Short-term PRMT5i",
                                                 "P3" = "Resistant PRMT5i",
                                                 "U1" = "Vehicle",
                                                 "U2" = "Vehicle",
                                                 "U3" = "Vehicle")
# Differential abundance  http://bioconductor.org/books/3.13/OSCA.multisample/differential-abundance.html

# long-term vs short-term AND untreated together
bb_cluster_representation2(
  obj = cds_human_pass_sf,
  sample_var = "sample",
  cluster_var = "leiden",
  comparison_var = "treatment_1",
  comparison_levels = c("Short-term PRMT5i/Vehicle", "Resistant PRMT5i"),
  sig_val = "FDR",
  return_val = "plot" # use "data" to get the table instead of the plot
)

# PRMT5i reistant vs untreated ONLY.
# Short-term treated sample filtered out
bb_cluster_representation2(
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

# resistant vs short_term treated ONLY.
# untreated samples filtered out
bb_cluster_representation2(
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

