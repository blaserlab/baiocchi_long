bb_rowmeta(cds_human_pass_sf) |> glimpse()
bb_gene_umap(cds_human_pass_sf,
             bb_rowmeta(cds_human_pass_sf) |>
               filter(mtor) |>
               mutate(label = "MTOR Pathway") |>
               select(feature_id, label), max_expr_val = 4)

bb_gene_umap(cds_human_pass_sf,
             bb_rowmeta(cds_human_pass_sf) |>
               filter(IGF1_signaling) |>
               mutate(label = "IGF1 Signaling") |>
               select(feature_id, label), max_expr_val = 3)

bb_gene_umap(cds_human_pass_sf,
             bb_rowmeta(cds_human_pass_sf) |>
               filter(eif4_p70S6k) |>
               mutate(label = "P70S6K Signaling") |>
               select(feature_id, label))

bb_gene_umap(cds_human_pass_sf,
             bb_rowmeta(cds_human_pass_sf) |>
               filter(ERK_MAPK_signaling) |>
               mutate(label = "ERK_MAPK") |>
               select(feature_id, label), max_expr_val = 3
             )


bb_gene_umap(cds_human_pass_sf,
             bb_rowmeta(cds_human_pass_sf) |>
               filter(pi3k_akt) |>
               mutate(label = "PI3K_AKT") |>
               select(feature_id, label), max_expr_val = 6
)
