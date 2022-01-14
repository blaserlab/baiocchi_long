# graphical parameters####
theme_font_size <- 12
theme_set(theme_cowplot(font_size = theme_font_size))


# show_col(pal_npg("nrc")(10))
experimental_group_palette <- brewer.pal(n = 8, name = "Set1")
# names(experimental_group_palette) <- analysis_configs$sample

jitter_alpha_fill <- 0.2
jitter_shape <- 21
jitter_size <- 2
jitter_stroke <- 0.5
jitter_width <- 0.2
jitter_alpha_color <- 1
jitter_height <- 0.2

summarybox_color <- "black"
summarybox_size <- 0.5
summarybox_width <- 0.3
summarybox_alpha <- 0.3
summarybox_geom <- "crossbar"

# 3 color heatmap
heatmap_3_colors <- c("#313695","white","#A50026")

# unmask
filter <- dplyr::filter
select <- dplyr::select
rename <- dplyr::rename
count <- dplyr::count

# output directories
figs_out <- "~/network/X/Labs/Blaser/collaborators/baiocchi_long_manuscript/figs/R_output"
tables_out <- "~/network/X/Labs/Blaser/collaborators/baiocchi_long_manuscript/tables"
