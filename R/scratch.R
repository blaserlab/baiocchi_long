# to make a new cell metadata column
colData(cds_human_pass_counts)$specimen_new <- recode(colData(cds_human_pass_counts)$specimen,
                                                      "P5" = "new_P5",
                                                      "" = "",
                                                      "" = "",
                                                      "" = "")
# make a faceted plot
bb_var_umap(cds_human_pass_counts, "specimen", value_to_highlight = c("P5", "P10", "P53", "P61")) + facet_wrap(facets = vars(value))

# make the overlay plot with custom color
bb_var_umap(cds_human_pass_counts, "specimen", value_to_highlight = c("P5"), palette = "blue")


colData(cds_human_pass_counts_5_10)

test_cds <- bb_align(cds_human_pass_counts_5_10, align_by = "specimen", n_cores = 24)
test_cds <- reduce_dimension(test_cds, cores = 39)
colData(test_cds)

bb_var_umap(test_cds, "specimen")
bb_var_umap(test_cds, "specimen", alt_dim_x = "prealignment_dim1", alt_dim_y = "prealignment_dim2")
identical(
  counts(cds_human_pass_counts_5_10),
  counts(test_cds)
)
bb_var_umap(test_cds, "partition")

test_cds <- bb_triplecluster(test_cds, n_top_markers = 50, outfile = "test.csv", n_cores = 39)

bb_var_umap(test_cds, "partition")
bb_var_umap(test_cds, "leiden")
bb_var_umap(test_cds, "louvain")
