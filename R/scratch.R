colData(cds_human_pass_sf)$erc_weeks <- factor(colData(cds_human_pass_sf)$erc_weeks, levels = c("0", "2", "8", "13"))
bb_var_umap(cds_human_pass_sf[,colData(cds_human_pass_sf)$tissue %in% c("Spleen")],
            "sample",
            facet_by = c("treatment","erc_weeks"),
            rows = vars(treatment),
            cols = vars(erc_weeks),
            foreground_alpha = 0.4,
            plot_title = "ERC x Treatment",
            palette = experimental_group_palette[str_detect(names(experimental_group_palette), "BM", negate = T)],
            legend_pos = "right")
