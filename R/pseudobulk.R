# not the preferred way
treatment_pseudosamples <- bb_cellmeta(cds_human_pass_sf) %>%
  group_by(sample, treatment, date_collected) %>%
  summarise()
treatment_pseudosamples

treatment_pseudobulk <-
  bb_pseudobulk_mf(
    cds = cds_human_pass_sf,
    pseudosample_table = treatment_pseudosamples,
    design_formula = "~date_collected+treatment",
    result_recipe = c("treatment", "PRMT5i", "Untreated")
  )


# the recommended way
partition_pseudosamples <- bb_cellmeta(cds_human_pass_sf) %>%
  group_by(sample, date_collected, partition, passage) %>%
  summarise()
partition_pseudosamples

partition_pseudobulk <-
  bb_pseudobulk_mf(
    cds = cds_human_pass_sf,
    pseudosample_table = partition_pseudosamples,
    design_formula = "~date_collected+passage+partition",
    result_recipe = c("partition", "1", "2")
  )
partition_pseudobulk$Result %>%
  filter(padj<0.01) %>%
  filter(log2FoldChange>0) %>%
  arrange(padj)

bb_gene_umap(cds_human_pass_sf,
             gene_or_genes = partition_pseudobulk$Result %>%
               filter(log2FoldChange > 0 & padj<0.05) %>%
               select(id) %>%
               mutate(gene_group = "up in PRMT5i"))

msigdb_geneset_metadata


partition_gsea <- fgsea::fgseaMultilevel(
  stats = deframe(
    partition_pseudobulk$Result %>% select(gene_short_name, log2FoldChange)
  ),
  pathways = bb_extract_msig(
    filter_list = list("ORGANISM" = "Homo sapiens"),
    return_form = "name_list"
  )
)
partition_gsea %>% filter(padj <0.05) %>% arrange(NES) %>% View()
