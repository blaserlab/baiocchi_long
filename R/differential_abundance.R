colData(cds_human_pass_sf)$treatment_1 <- recode(colData(cds_human_pass_sf)$sample,
                                                 "P1" = "long-term",
                                                 "P2" = "short-term/untreated",
                                                 "P3" = "long-term",
                                                 "U1" = "short-term/untreated",
                                                 "U2" = "short-term/untreated",
                                                 "U3" = "short-term/untreated")
bb_var_umap(cds_human_pass_sf, "treatment_1", facet_by = "value")

# Differential abundance  http://bioconductor.org/books/3.13/OSCA.multisample/differential-abundance.html
leiden_counts <- bb_cellmeta(cds_human_pass_sf) |>
  filter(treatment == "PRMT5i") |>
  count(sample, leiden) |>
  pivot_wider(names_from = "sample", values_from = "n", values_fill = 0) |>
  bb_tbl_to_matrix()
leiden_counts

dge_list <- edgeR::DGEList(counts = leiden_counts,
                           samples = bb_cellmeta(cds_human_pass_sf) |>
                             filter(treatment == "PRMT5i") |>
                             group_by(sample,
                                      treatment_1) |>
                             summarise())

design <- model.matrix(~factor(treatment_1), dge_list$samples)
dge_list <- edgeR::estimateDisp(dge_list, design, trend="none")
fit <- edgeR::glmQLFit(dge_list, design, robust=TRUE, abundance.trend=FALSE)
res <- edgeR::glmQLFTest(fit, coef=ncol(design))
res_top_tags <- edgeR::topTags(res, n = Inf)
cluster_enrichment_barchart <- as_tibble(res_top_tags@.Data[[1]], rownames = "cluster") |>
  mutate(logFC = -logFC) |>
  mutate(enriched = ifelse(logFC < 0, "short-term/untreated", "long-term")) |>
  mutate(sig = case_when(PValue < 0.05 & PValue >= 0.01 ~ "*",
                         PValue < 0.01 & PValue >= 0.001 ~ "**",
                         PValue < 0.001 & PValue >= 0.0001 ~ "***",
                         PValue < 0.0001 ~ "****",
                         PValue >= 0.05 ~ ""
                         )) |>
  mutate(cluster = fct_reorder(cluster, logFC)) |>
  mutate(texty = ifelse(
    logFC > 0,
    logFC,
    0
  )) |>
  ggplot(aes(x = cluster,
             y = logFC,
             fill = logFC,
             label = sig)) +
  scale_fill_viridis_c() +
  geom_col() +
  labs(x = NULL,
       # y = "Log<sub>2</sub>Fold Enrichment:<br>Long-Term Treated/Short-Term and Untreated",
       y = "Log<sub>2</sub>Fold Enrichment:<br>Long-Term Treated/Short-Term Treated",
       color = NULL,
       fill = NULL) +
  theme(axis.title.y.left = ggtext::element_markdown()) +
  geom_text(aes(y = texty), color = "black", nudge_y = 0.05) +
  theme(legend.position = "none")
cluster_enrichment_barchart
