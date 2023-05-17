library("data.table")
library("tidyverse")
library("ggplot2")
library("glue")
library("ggpubr")

options(scipen=10000)

counts = "../data/count_tables/20230510_scTIP_read_sums.tsv"
meta = "../data/scTIP_seq_wells_230305.txt"

ct = fread(glue("{counts}"))
ct = ct %>% mutate(pos = unname(sapply(ct$well, function(x) strsplit(x, split = "_")[[1]][2])))
ct = ct %>% mutate(sample_id = as.numeric(unname(sapply(ct$well, function(x) strsplit(x, split = "_")[[1]][3])))) %>% 
  mutate(sample = case_when((sample_id >= 1 & sample_id <= 96) ~ "EpiLC_6h",
            (sample_id > 96 & sample_id <= 192) ~ "EpiLC_12h",
            (sample_id > 192 & sample_id <= 288) ~ "EpiLC_24h",
            (sample_id > 288 & sample_id <= 384) ~ "EpiLC_48h", 
            TRUE ~ "-"))

summary(ct$overal_read_count)

total_rc_histo = ggplot(ct, aes(x = overal_read_count)) +
  geom_histogram(position = "identity", fill = "#fc9272") +
  labs(
    title = "",
    x = "",
    y = "total read counts"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black")
  )
total_rc_histo

sample_histo = ggplot(ct, aes(x = overal_read_count, fill = sample)) +
  geom_histogram(position = "identity", alpha=0.5) +
  scale_fill_brewer(palette = "Set3") +
  labs(
    title = "",
    x = "",
    y = "total read counts"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black")
  )
sample_histo

ggsave(
  "../results/20230510_total_rc_samples_histogram.png",
  plot = last_plot(),
  height = 6,
  width = 6,
  dpi = 500
)

ggsave(
  "../results/20230510_total_rc_samples_histogram.pdf",
  plot = last_plot(),
  height = 6,
  width = 6
)

order = factor(ct$sample, levels = c("EpiLC_6h", "EpiLC_12h", "EpiLC_24h", "EpiLC_48h"))
sample_box = ggplot(ct, aes(x = order, y = overal_read_count, fill = sample)) +
  geom_boxplot(color = "black") +
  scale_fill_brewer(palette = "Set3") +
  labs(
    title = "",
    x = "",
    y = "total read counts",
    fill = ""
  ) +
  guides(fill = "none") +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 20, color = "black", angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 20, color = "black")
  )
sample_box

ggsave(
  "../results/20230510_total_read_count_box.png",
  plot = last_plot(),
  height = 6,
  width = 6,
  dpi = 500
)

ggsave(
  "../results/20230510_total_read_count_box.pdf",
  plot = last_plot(),
  height = 6,
  width = 6
)


