library("data.table")
library("tidyverse")
library("ggplot")
library("glue")
library("ggpubr")

counts = "../data/count_tables/2023apr_scTIP_read_sums.tsv"
meta = "../data/scTIP_seq_wells_230305.txt"

ct = fread(glue("{counts}"))
ct = ct %>% mutate(pos = unname(sapply(ct$well, function(x) strsplit(x, split = "_")[[1]][2])))

meta = fread(meta)
meta = meta %>% separate(well, "_", into = c("sample", "pos")) %>% mutate(sample = ifelse(sample == "H33", "H3.3", sample))

ct_meta = ct %>% inner_join(., meta, by = "pos")

ggplot(ct_meta, aes(x = sample, y = overal_read_count, fill = sample)) +
  geom_boxplot(color = "black") +
  scale_fill_brewer(palette = "Reds") +
  ylim(0, 100) +
  labs(
    title = "",
    x = "",
    y = "total read count / well",
    fill = ""
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black")
  ) +
  stat_compare_means(label.y = 75, label.x = 1.25, size = 6) 
  # stat_compare_means(label = "p.signif", method = "t.test",
  #                    ref.group = ".all.", label.y = 50)


ggsave(
  "../results/total_read_counts_per_well.png",
  plot = last_plot(),
  height = 7,
  width = 7,
  dpi = 500
)

ggsave(
  "../results/total_read_counts_per_well.pdf",
  plot = last_plot(),
  height = 7,
  width = 7
)

above_100 = ct_meta %>% mutate(above_100 = ifelse(100 < overal_read_count, "> 100", "< 100")) %>% 
  group_by(above_100) %>% count() %>% 
  ggplot(data = ., aes(x = above_100, y = n)) +
  geom_bar(stat = "identity", fill = "steelblue")+
  geom_text(aes(label = n), vjust=1.6, color="white", size=7.5) +
  labs(
    title = "",
    x = "read count",
    y = "# of cells",
    fill = ""
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black")
  )
above_100

ggsave(
  "../results/above_100_bar.png",
  plot = above_100,
  height = 7,
  width = 7,
  dpi = 500
)

ggsave(
  "../results/above_100_bar.pdf",
  plot = above_100,
  height = 7,
  width = 7
)
