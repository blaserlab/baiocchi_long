colData(cds_human_pass_sf)$treatment_1 <- recode(colData(cds_human_pass_sf)$sample,
                                                 "P1" = "long-term",
                                                 "P2" = "short-term/untreated",
                                                 "P3" = "long-term",
                                                 "U1" = "short-term/untreated",
                                                 "U2" = "short-term/untreated",
                                                 "U3" = "short-term/untreated")
bb_var_umap(cds_human_pass_sf, "treatment_1", facet_by = "value")

# Differential abundance  http://bioconductor.org/books/3.13/OSCA.multisample/differential-abundance.html
#params
obj <- cds_human_pass_sf
cluster_var <- "leiden"
sample_var <- "sample"
comparison_var <- "treatment_1"
comparison_levels <- c("short-term/untreated", "long-term") # clusters enriched in the first term will be negative and the second term will be positive

# get counts by cluster
counts <- bb_cellmeta(obj) |>
  dplyr::count(!!sym(sample_var), !!sym(cluster_var)) |>
  tidyr::pivot_wider(names_from = sample_var,
              values_from = "n",
              values_fill = 0) |>
  bb_tbl_to_matrix()

dge_list <- edgeR::DGEList(counts = counts,
                           samples = bb_cellmeta(obj) |>
                             group_by(!!sym(sample_var),
                                      !!sym(comparison_var)) |>
                             summarise())

dge_list$samples[[comparison_var]] <- factor(dge_list$samples[[comparison_var]], levels = comparison_levels)
design <- model.matrix(~dge_list$samples[[comparison_var]], dge_list$samples,)
dge_list <- edgeR::estimateDisp(dge_list, design, trend="none")
fit <- edgeR::glmQLFit(dge_list, design, robust=TRUE, abundance.trend=FALSE)
res <- edgeR::glmQLFTest(fit, coef=ncol(design))
res_top_tags <- edgeR::topTags(res, n = Inf)
data <- as_tibble(res_top_tags@.Data[[1]], rownames = "cluster") |>
  mutate(enriched = ifelse(logFC < 0, comparison_levels[1], comparison_levels[2])) |>
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
  ))
plot <-
  ggplot(data, aes(x = cluster,
             y = logFC,
             fill = enriched,
             label = sig)) +
  geom_col(color = "black") +
  scale_fill_manual(values = experimental_group_palette_1) +
  # scale_fill_viridis_c() +
  # scale_fill_gradient2(low = "red3", mid = "white", high = "blue4", midpoint = 0) +
  labs(x = cluster_var,
       y = "Log<sub>2</sub>Fold Enrichment",
       fill = "Enriched Group") +
  theme(axis.title.y.left = ggtext::element_markdown()) +
  geom_text(aes(y = texty), color = "black", nudge_y = 0.05) +
  theme(legend.position = "right")
