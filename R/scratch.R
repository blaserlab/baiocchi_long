system.file("data-raw/load_process.R", package = "baiocchi.long.datapkg")

# navigate to renv/library/R-4.2/x86_64-pc-linux-gnu/baiocchi.long.datapkg/data-raw/load_process.R

"~/network/X/Labs/Blaser/single_cell/baiocchi_long_20211209/analysis_config.csv"
# make the directory
dir.create(file.path(figs_out, "mac_geneset_figs"))

# produce the plots
walk2(
  @@ -8,9 +12,15 @@ walk2(
    "mtor",
    "pi3k_akt",
    "pi3k_bcells",
    "eif4_p70S6k"
    "eif4_p70S6k",
    "ERK_MAPK_signaling",
    "IGF1_signaling",
    "insulin_secretion_signaling",
    "LPS_MAPK_signaling",
    "UVA_MAPK_signaling"

  ),
  .y = c(Inf, 4, 6, 10, Inf),
  .y = c(Inf, 4, 6, 10, Inf, 3.5, Inf, 2.5, 4, 3),
  .f = \(x, y, dat = cds_human_pass_sf, figdir = figs_out) {
    p <-
      bb_gene_umap(
        @@ -34,12 +44,19 @@ walk2(

          # print out the gene sets used
          bb_rowmeta(cds_human_pass_sf) |>
            select(gene_short_name,
                   ER_signaling,
                   mtor,
                   pi3k_akt,
                   pi3k_bcells,
                   eif4_p70S6k) |>
            select(
              gene_short_name,
              ER_signaling,
              mtor,
              pi3k_akt,
              pi3k_bcells,
              eif4_p70S6k,
              ERK_MAPK_signaling,
              IGF1_signaling,
              insulin_secretion_signaling,
              LPS_MAPK_signaling,
              UVA_MAPK_signaling
            ) |>
            pivot_longer(-gene_short_name, names_to = "pathway") |>
            filter(value) |>
            select(gene_short_name, pathway) |>
