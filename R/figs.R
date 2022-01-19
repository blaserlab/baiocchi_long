# UMAPs
bb_var_umap(cds_human_pass_sf, "sample")
bb_var_umap(cds_human_pass_sf, "orig_id", facet_by = "value")

treatment_density_umap <- bb_var_umap(cds_human_pass_sf, "density", facet_by = "treatment")
save_plot(treatment_density_umap, filename = str_glue("{figs_out}/treatment_density_umap.png"), base_width = 6.5, base_height = 3.5)

bb_var_umap(cds_human_pass_sf, "clonotype")
bb_var_umap(cds_human_pass_sf, "partition")
bb_var_umap(cds_human_pass_sf, "treatment")
bb_var_umap(cds_human_pass_sf, "leiden")
bb_var_umap(cds_human_pass_sf, "leiden",
            facet_by = c("treatment", "orig_id"), rows = vars(treatment), cols = vars(orig_id))

bb_var_umap(cds_human_pass_sf, "leiden") +
  facet_grid(cols = vars(orig_id), rows = vars(treatment))

# get mackenzies genes for plotting

mackenzie_genesets <-
  c(
    "BIOCARTA_ERK_PATHWAY",
    "REACTOME_FCERI_MEDIATED_NF_KB_ACTIVATION",
    "HAMAI_APOPTOSIS_VIA_TRAIL_UP",
    "FISCHER_G1_S_CELL_CYCLE",
    "GO_SMAD_PROTEIN_SIGNAL_TRANSDUCTION",
    "DACOSTA_UV_RESPONSE_VIA_ERCC3_COMMON_DN",
    "DACOSTA_UV_RESPONSE_VIA_ERCC3_TTD_DN",
    "VANASSE_BCL2_TARGETS_UP",
    "FISCHER_DREAM_TARGETS",
    "NAGASHIMA_EGF_SIGNALING_UP",
    "BIOCARTA_BCELLSURVIVAL_PATHWAY",
    "BIOCARTA_GPCR_PATHWAY",
    "BIOCARTA_ARENRF2_PATHWAY"
  )


# blaseRdata::msigdb_geneset_metadata
mackenzie_geneset_data <- bb_extract_msig(filter_list = list("STANDARD_NAME" = mackenzie_genesets), return_form = "tibble") %>%
  filter(!is.na(gene_short_name)) %>%
  select(-id) %>%
  left_join(bb_rowmeta(cds_human_pass_sf)) %>%
  select(id, gene_set) %>%
  filter(!is.na(id)) %>%
  group_by(gene_set, id) %>%
  summarise() %>%
  relocate(id)
# mackenzie_geneset_data <- bind_rows(mackenzie_geneset_data, non_msig_data)



geneset_plots <- map(.x = mackenzie_geneset_data %>%
                       pull(gene_set) %>%
                       unique(),
    .f = function(x,
                  cds = cds_human_pass_sf,
                  data = mackenzie_geneset_data) {
      genes_to_plot <- data %>%
        filter(gene_set == x)
      p <- bb_gene_umap(cds, gene_or_genes = genes_to_plot, plot_title = x) +
        facet_wrap(facets = vars(treatment))
    })


names(geneset_plots) <- mackenzie_geneset_data %>% pull(gene_set) %>% unique()

geneset_plots$BIOCARTA_ARENRF2_PATHWAY
geneset_plots$BIOCARTA_BCELLSURVIVAL_PATHWAY
geneset_plots$BIOCARTA_ERK_PATHWAY
geneset_plots$BIOCARTA_GPCR_PATHWAY
geneset_plots$DACOSTA_UV_RESPONSE_VIA_ERCC3_COMMON_DN
geneset_plots$DACOSTA_UV_RESPONSE_VIA_ERCC3_TTD_DN
geneset_plots$FISCHER_DREAM_TARGETS
geneset_plots$FISCHER_G1_S_CELL_CYCLE
geneset_plots$HAMAI_APOPTOSIS_VIA_TRAIL_UP
geneset_plots$NAGASHIMA_EGF_SIGNALING_UP
geneset_plots$REACTOME_FCERI_MEDIATED_NF_KB_ACTIVATION
geneset_plots$VANASSE_BCL2_TARGETS_UP

bb_gene_umap(cds_human_pass_sf, gene_or_genes = bb_rowmeta(cds_human_pass_sf) %>% select(id, module))

bb_gene_umap(cds_human_pass_sf, gene_or_genes = "MTOR")
bb_gene_umap(cds_human_pass_sf, gene_or_genes = "MTOR") + scale_color_viridis_c()
bb_gene_umap(cds_human_pass_sf, gene_or_genes = "MIR34AHG") + scale_color_viridis_c()

# cluster representataion
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
  cds = cds_human_pass_sf[,colData(cds_human_pass_sf)$treatment_time %in% c("short_term_treated", "Untreated")],
  cluster_var = "partition",
  class_var = "treatment_time",
  experimental_class = "short_term_treated",
  control_class = "Untreated",
  return_value = "table"
)


bb_cluster_representation(
  cds = cds_human_pass_sf[,colData(cds_human_pass_sf)$treatment_time %in% c("long_term_treated", "Untreated")],
  cluster_var = "partition",
  class_var = "treatment_time",
  experimental_class = "long_term_treated",
  control_class = "Untreated",
  return_value = "table"
)
bb_cluster_representation(
  cds = cds_human_pass_sf,
  cluster_var = "leiden_3_binary",
  class_var = "treatment",
  experimental_class = "PRMT5i",
  control_class = "Untreated",
  return_value = "plot"
)


bb_cellmeta(cds_human_pass_sf) %>%
  group_by(sample, leiden_3_binary) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = "leiden_3_binary", values_from = "n") %>%
  mutate(ratio = leiden_3/not_leiden_3) %>%
  mutate(log2ratio = log2(ratio))
