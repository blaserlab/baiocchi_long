coldata_with_counts <- left_join(
  colData(cds_human) %>%
    as_tibble(rownames = "column"),
  cds_human_counts
)

coldata_with_counts %>%
  group_by(partition) %>%
  summarise(mean_sf = mean(Size_Factor), mean_counts = mean(sum_counts))

tibble(sum_counts = density(coldata_with_counts$sum_counts)$x,
density = density(coldata_with_counts$sum_counts)$y) %>%
ggplot(coldata_with_counts, mapping = aes(x = sum_counts, y = density)) +
  geom_line() +
  ggpmisc::stat_valleys(colour = "red") +
  ggpmisc::stat_valleys(geom = "text", hjust = -0.2, vjust = 0.1, angle = 45)

tibble(Size_Factor = density(coldata_with_counts$Size_Factor)$x,
density = density(coldata_with_counts$Size_Factor)$y) %>%
ggplot(coldata_with_counts, mapping = aes(x = Size_Factor, y = density)) +
  geom_line() +
  ggpmisc::stat_valleys(colour = "red") +
  ggpmisc::stat_valleys(geom = "text", hjust = -0.2, vjust = 0.1, angle = 45)
