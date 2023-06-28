library("data.table")
library("GenomicRanges")
library("tidyverse")
library("wigglescout")
library("ggrastr")
library("ggpubr")

peaks_6h = fread("../results/SEACR/EpiLC_6h.stringent_thr0.2.bed")
peaks_12h = fread("../results/SEACR/EpiLC_12h.stringent_thr0.2.bed")
peaks_24h = fread("../results/SEACR/EpiLC_24h.stringent_thr0.2.bed")
peaks_48h = fread("../results/SEACR/EpiLC_48h.stringent_thr0.2.bed")

# Cruz-Molina active enhancers
cm = fread("../data/bed/ESC_Enhancer_CruzMolina.active_mm10.bed")

gr = function(peaks, sample) {
  peaks$sample = sample
  peaks = GRanges(
    seqnames = peaks$V1,
    ranges = IRanges(
      start = peaks$V2,
      end = peaks$V3,
      score = peaks$V5,
      names = peaks$sample,
    )
  )
  
  return(peaks)
}

peaks_6h_gr = gr(peaks_6h, sample = "6h")
peaks_12h_gr = gr(peaks_12h, sample = "12h")
peaks_24h_gr = gr(peaks_24h, sample = "24h")
peaks_48h_gr = gr(peaks_48h, sample = "48h")

peaks = GRangesList("6h" = peaks_6h_gr, "12h" = peaks_12h_gr, "24h" = peaks_24h_gr, "48h" = peaks_48h_gr)

bigwigs = c("../data/20230510_EpiLC/EpiLC_6h/EpiLC_6h_pseudobulk_RPGC.bigwig",
            "../data/20230510_EpiLC/EpiLC_12h/EpiLC_12h_pseudobulk_RPGC.bigwig",
            "../data/20230510_EpiLC/EpiLC_24h/EpiLC_24h_pseudobulk_RPGC.bigwig",
            "../data/20230510_EpiLC/EpiLC_48h/EpiLC_48h_pseudobulk_RPGC.bigwig")

read_counts = bw_loci(bigwigs, Reduce(subsetByOverlaps, peaks))
read_counts = as.data.frame(read_counts)

mat = read_counts %>% dplyr::select(starts_with("Epi")) %>% as.matrix

monotonic_incr = function(signal) {
  return(all(signal == cummax(signal)))
}

monotonic_decreas = function(signal) {
  return(all(signal == cummin(signal)))
}


increasing_H33 = read_counts[which(apply(mat, MARGIN = 1, monotonic_incr)),]
increasing_H33 = increasing_H33 %>% mutate(diff = EpiLC_48h_pseudobulk_RPGC - EpiLC_6h_pseudobulk_RPGC)
increasing_H33_long = increasing_H33 %>% pivot_longer(., cols = "EpiLC_6h_pseudobulk_RPGC":"EpiLC_48h_pseudobulk_RPGC",
                                                      values_to = "read_count", names_to = "sample")

increasing_H33_long = increasing_H33_long %>% 
  mutate(id = paste(seqnames, start, end, sep = "_"))

ids = unique(increasing_H33_long$id)
increasing_H33_long = lapply(ids, function(x) {increasing_H33_long %>% filter(id == x) %>% mutate(order = seq(1, 4))})
increasing_H33_long = bind_rows(increasing_H33_long)

smooth1 = ggplot() + geom_smooth(
  data = increasing_H33_long,
  aes(x = order, y = read_count),
  se = TRUE,
  color = "#9ecae1",
  method = "lm",
  formula = "y ~ x",
  linewidth = 0.5
) +
  labs(
    title = "monotonic increasing H3.3 signal",
    x = "",
    y = "read count"
  ) +
  ylim(0, 3) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black")
  ) 

increasing_H33_long = increasing_H33_long %>%
  mutate(across('sample', str_replace, 'EpiLC_6h_pseudobulk_RPGC', '6h')) %>%
  mutate(across('sample', str_replace, 'EpiLC_12h_pseudobulk_RPGC', '12h')) %>%
  mutate(across('sample', str_replace, 'EpiLC_24h_pseudobulk_RPGC', '24h')) %>%
  mutate(across('sample', str_replace, 'EpiLC_48h_pseudobulk_RPGC', '48h'))

sample_order = factor(increasing_H33_long$sample, levels = c("6h", "12h", "24h", "48h"))

boxplot1 = ggplot() + geom_boxplot(
  data = decreasing_H33_long,
  aes(x = sample_order, y = read_count),
  color = "#9ecae1"
) +
  labs(
    title = "",
    x = "time point",
    y = "read count"
  ) +
  ylim(0, 50) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black")
  ) 


decreasing_H33 = read_counts[which(apply(mat, MARGIN = 1, monotonic_decreas)),]
decreasing_H33 = decreasing_H33 %>% mutate(diff = EpiLC_6h_pseudobulk_RPGC - EpiLC_48h_pseudobulk_RPGC)
decreasing_H33_long = decreasing_H33 %>% pivot_longer(., cols = "EpiLC_6h_pseudobulk_RPGC":"EpiLC_48h_pseudobulk_RPGC",
                                                      values_to = "read_count", names_to = "sample")

decreasing_H33_long = decreasing_H33_long %>% 
  mutate(id = paste(seqnames, start, end, sep = "_"))

ids = unique(decreasing_H33_long$id)
decreasing_H33_long = lapply(ids, function(x) {decreasing_H33_long %>% filter(id == x) %>% mutate(order = seq(1, 4))})
decreasing_H33_long = bind_rows(decreasing_H33_long)

smooth2 = ggplot() + geom_smooth(
  data = decreasing_H33_long,
  aes(x = order, y = read_count),
  se = TRUE,
  color = "#fc9272",
  method = "lm",
  formula = "y ~ x",
  linewidth = 0.5
) +
  labs(
    title = "monotonic decreasing H3.3 signal",
    x = "",
    y = "read count"
  ) +
  ylim(0, 3) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black")
  ) 

decreasing_H33_long = decreasing_H33_long %>%
  mutate(across('sample', str_replace, 'EpiLC_6h_pseudobulk_RPGC', '6h')) %>%
  mutate(across('sample', str_replace, 'EpiLC_12h_pseudobulk_RPGC', '12h')) %>%
  mutate(across('sample', str_replace, 'EpiLC_24h_pseudobulk_RPGC', '24h')) %>%
  mutate(across('sample', str_replace, 'EpiLC_48h_pseudobulk_RPGC', '48h'))
  
sample_order = factor(decreasing_H33_long$sample, levels = c("6h", "12h", "24h", "48h"))

boxplot2 = ggplot() + geom_boxplot(
  data = decreasing_H33_long,
  aes(x = sample_order, y = read_count),
  color = "#fc9272"
) +
  labs(
    title = "",
    x = "time point",
    y = "read count"
  ) +
  ylim(0, 50) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black")
  ) 

gr1 = read_counts
gr1$type = "H3.3"
gr1 = GRanges(
  seqnames = gr1$seqnames,
  ranges = IRanges(
    start = gr1$start,
    end = gr1$end,
    names = gr1$type,
  )
)
cm$type = "active_enh"
gr2 = GRanges(
  seqnames = cm$V1,
  ranges = IRanges(
    start = cm$V2,
    end = cm$V3,
    names = cm$type,
  )
)

ol = findOverlaps(gr1, gr2, type = "any", ignore.strand = FALSE)
h33_enh = as.data.frame(read_counts[queryHits(ol),])
h33_enh = h33_enh %>% 
  mutate(row_id = queryHits(ol)) %>% 
  mutate(enhancer_status = "active enhancer") %>% 
  dplyr::select(enhancer_status, row_id)
#read_counts = read_counts %>% rowid_to_column(., "row_id")
h33_enh = read_counts %>% left_join(., h33_enh, by = "row_id") %>% 
  dplyr::mutate(enhancer_status = replace_na(enhancer_status, "non-enhancer"))
mat = h33_enh %>% dplyr::select(starts_with("Epi")) %>% as.matrix

h33_enh_incr = h33_enh %>% mutate(monotonicity = apply(mat, MARGIN = 1, monotonic_incr)) %>% 
  mutate(monotonicity = ifelse(monotonicity == TRUE, "increase", NA_character_)) %>% 
  dplyr::filter(monotonicity == "increase")

h33_enh_incr_long = h33_enh_incr %>% pivot_longer(., cols = "EpiLC_6h_pseudobulk_RPGC":"EpiLC_48h_pseudobulk_RPGC",
                                                      values_to = "read_count", names_to = "sample")

h33_enh_incr_long = h33_enh_incr_long %>% 
  mutate(id = paste(seqnames, start, end, sep = "_"))

blacklist = h33_enh_incr_long %>% group_by(id) %>% count(id) %>% dplyr::filter(n != 4) %>% pull(id)
h33_enh_incr_long = h33_enh_incr_long %>% dplyr::filter(!id %in% blacklist)
ids = unique(h33_enh_incr_long$id)
h33_enh_incr_long = lapply(ids, function(x) {h33_enh_incr_long %>% filter(id == x) %>% mutate(order = seq(1, 4))})
h33_enh_incr_long = bind_rows(h33_enh_incr_long)

smooth3 = ggplot() + geom_smooth(
  data = h33_enh_incr_long,
  aes(x = order, y = read_count, fill = enhancer_status),
  se = TRUE,
  method = "lm",
  formula = "y ~ x",
  linewidth = 0.5
) +
  scale_fill_manual(values = c("#fc9272", "#9ecae1")) +
  labs(
    title = "Increasing signal over Cruz-Molina enhancers",
    x = "",
    y = "read count",
    fill = "status"
  ) +
  ylim(0, 3) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black")
  ) 
smooth3

h33_enh_decr = h33_enh %>% mutate(monotonicity = apply(mat, MARGIN = 1, monotonic_decreas)) %>% 
  mutate(monotonicity = ifelse(monotonicity == TRUE, "decrease", NA_character_)) %>% 
  dplyr::filter(monotonicity == "decrease")

h33_enh_decr_long = h33_enh_decr %>% pivot_longer(., cols = "EpiLC_6h_pseudobulk_RPGC":"EpiLC_48h_pseudobulk_RPGC",
                                                  values_to = "read_count", names_to = "sample")

h33_enh_decr_long = h33_enh_decr_long %>% 
  mutate(id = paste(seqnames, start, end, sep = "_"))

blacklist = h33_enh_decr_long %>% group_by(id) %>% count(id) %>% dplyr::filter(n != 4) %>% pull(id)
h33_enh_decr_long = h33_enh_decr_long %>% dplyr::filter(!id %in% blacklist)
ids = unique(h33_enh_decr_long$id)
h33_enh_decr_long = lapply(ids, function(x) {h33_enh_decr_long %>% filter(id == x) %>% mutate(order = seq(1, 4))})
h33_enh_decr_long = bind_rows(h33_enh_decr_long)

smooth4 = ggplot() + geom_smooth(
  data = h33_enh_decr_long,
  aes(x = order, y = read_count, fill = enhancer_status),
  se = TRUE,
  method = "lm",
  formula = "y ~ x",
  linewidth = 0.5
) +
  scale_fill_manual(values = c("#fc9272", "#9ecae1")) +
  labs(
    title = "Decreasing signal over Cruz-Molina enhancers",
    x = "",
    y = "read count",
    fill = "status"
  ) +
  ylim(0, 3) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black")
  ) 
smooth4

ggarrange(smooth1, smooth2, boxplot1, boxplot2, smooth3, smooth4, ncol = 2, nrow = 3)

ggsave(
  "../results/EpiLC_H3.3_monotonic_changes.png",
  plot = last_plot(),
  height = 8,
  width = 10,
  dpi = 500
)

ggsave(
  "../results/EpiLC_H3.3_monotonic_changes.pdf",
  plot = last_plot(),
  height = 8,
  width = 10
)
