if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse",
               "data.table",
               "wigglescout",
               "ggpubr",
               "glue",
               "RColorBrewer",
               "GenomicRanges",
               "magick",
               "ggrastr",
               "ComplexHeatmap",
               "circlize",
               "matrixStats") 


# output folder
result_folder = "../results/wigglescout/"

# EpiLC H3.3 sciTIP-Seq signals
epilc = list.files("../results/Seurat/cluster_bigwigs/H3.3_EpiLC_timepoints/",
                        pattern = "*bigwig",
                        full.names = TRUE)
epilc_peaks = list.files("../results/SEACR/",
                         pattern = "*bed",
                         full.names = TRUE)

# EpiLC ChIP-Seq (Yang et al.)
yang = list.files("../data/bigwig/Yang_2019/", pattern = "*.bw", full.names = TRUE)

# Cruz-Molina active enhancers (mESC)
cm = "../data/bed/ESC_Enhancer_CruzMolina.active_mm10.bed"
cm = fread(cm)
cm = cm %>% mutate(V2 = V2 - 1000, V3 = V3 + 1000)

cm$V6 = "Cruz-Molina_active_enh"
cm = GRanges(
  seqnames = cm$V1,
  ranges = IRanges(
    start = cm$V2,
    end = cm$V3,
    names = cm$V6
  )
)

cm_orig_size = "../data/bed/ESC_Enhancer_CruzMolina.active_mm10.bed"
cm_orig_size = fread(cm_orig_size)
cm_orig_size$V6 = "Cruz-Molina_active_enh"
cm_orig_size = GRanges(
  seqnames = cm_orig_size$V1,
  ranges = IRanges(
    start = cm_orig_size$V2,
    end = cm_orig_size$V3,
    names = cm_orig_size$V6
  )
)

## helper functions
# create scatter plot between two Minute-ChIP bigwigs
create_scatter = function(bigwig1 = epilc[1], bigwig2 = yang[1]) {
  name1 = strsplit(strsplit(bigwig1, "../results/Seurat/cluster_bigwigs/H3.3_EpiLC_timepoints/")[[1]][2],
                   "_pseudobulk_RPGC.bigwig")[[1]][1]
  name2 = strsplit(strsplit(bigwig2, "../data/bigwig/Yang_2019/")[[1]][2],
                   "_")[[1]][2]
  
  done = list.files(glue("{result_folder}"), full.names = FALSE)
  if (glue("{name1}-{name2}_sc.pdf") %in% done) {
    return(print(glue("{name1}-{name2} is done!")))
  } else {
    print(glue("Working on: {name1}-{name2}"))
  }
  
  peaks = epilc_peaks[grep(name1, epilc_peaks)]
  
  peaks = fread(peaks)
  peaks = GRanges(
    seqnames = peaks$V1,
    ranges = IRanges(
      start = peaks$V2,
      end = peaks$V3
    )
  )
  
  ol = findOverlaps(peaks, cm, type = "any", ignore.strand = TRUE)
  ol_w_cm = peaks[queryHits(ol)]
  
  p = plot_bw_loci_scatter(
    bigwig1,
    bigwig2,
    loci = peaks,
    highlight = ol_w_cm,
    highlight_colors = "#de2d26",
    highlight_label = "enh",
    remove_top = 0.0001,
    verbose = FALSE
  ) + 
    xlim(0, 20) +
    ylim(0, 20) +
    theme(
      text = element_text(size = 16),
      legend.text = element_text(size = 14),
      plot.title = element_text(size = 7, face = "bold"),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title.x = element_text(size = 14, color = "black"),
      axis.title.y = element_text(size = 14, color = "black")
    ) +
    guides(color = "none", fill = "none") +
    geom_abline(slope = 1,
                intercept = 0,
                color = "#de2d26",
                linetype = "dashed") +
    labs(title = glue("H3.3 {name1} - {name2}"),
         x = name1,
         y = name2) +
    stat_cor(method = "pearson", label.x = 3, label.y = 17, size = 2)
  # rasterize!
  p = rasterize(p, layers='Point', dpi = 150)
  # save
  ggsave(
    glue("{result_folder}{name1}-{name2}_sc.pdf"),
    plot = p,
    width = 3,
    height = 3
  )
  return(p)
}

# calculate correlations bitween read counts
# considering the top H3.3 peaks! (see thr)
correlation_analysis = function(bigwig1 = epilc[1],
                                bigwig2 = yang[1]) {
  name1 = strsplit(
    strsplit(
      bigwig1,
      "../results/Seurat/cluster_bigwigs/H3.3_EpiLC_timepoints/"
    )[[1]][2],
    "_pseudobulk_RPGC.bigwig"
  )[[1]][1]
  name2 = strsplit(strsplit(bigwig2, "../data/bigwig/Yang_2019/")[[1]][2],
                   "_")[[1]][2]
  
  peaks = epilc_peaks[grep(name1, epilc_peaks)]
  
  peaks = fread(peaks)
  peaks = GRanges(seqnames = peaks$V1,
                  ranges = IRanges(start = peaks$V2,
                                   end = peaks$V3))
  
  ol = findOverlaps(peaks, cm, type = "any", ignore.strand = TRUE)
  ol_w_cm = peaks[queryHits(ol)]
  
  print(glue("Working on: {name1}-{name2}"))
  read_densities = bw_loci(c(bigwig1, bigwig2),
                           loci = ol_w_cm)
  read_densities = as_tibble(read_densities)
  thr = read_densities %>% select(starts_with("Epi")) %>% arrange(., across(starts_with("Epi"), desc)) %>% pull(starts_with("Epi")) 
  thr = thr[50]
  read_densities = read_densities %>% select(-seqnames,-start,-end,-width,-strand) %>% 
    filter_at(1, all_vars(. > thr)) %>% 
    na.omit()
  epilc = read_densities %>% pull(starts_with("Epi"))
  yang = read_densities %>% pull(starts_with("GS"))
  pearson = round(cor(epilc, yang, method = "pearson"), 2)
  spearman = round(cor(epilc, yang, method = "spearman"), 2)
  
  comparison = glue("{name1}-{name2}")
  output = tibble(
    method = c("pearson", "spearman"),
    coefficient = c(pearson, spearman),
    comparison = rep(comparison, 2)
  )
  return(output)
}

top_h33_over_cm =  function(h33_bigwig = epilc[4], top = 50) {
  name = strsplit(
    strsplit(
      h33_bigwig,
      "../results/Seurat/cluster_bigwigs/H3.3_EpiLC_timepoints/"
    )[[1]][2],
    "_pseudobulk_RPGC.bigwig"
  )[[1]][1]
  
  peaks = epilc_peaks[grep(name, epilc_peaks)]
  
  peaks = fread(peaks)
  peaks = GRanges(seqnames = peaks$V1,
                  ranges = IRanges(start = peaks$V2,
                                   end = peaks$V3))
  
  ol = findOverlaps(peaks, cm, type = "any", ignore.strand = TRUE)
  ol_w_cm = peaks[queryHits(ol)]
  
  read_densities = bw_loci(h33_bigwig,
                           loci = ol_w_cm)
  read_densities = as_tibble(read_densities)
  thr = read_densities %>% pull(starts_with("Epi"))
  thr = thr[order(thr, decreasing = TRUE)][top]
  
  read_densities = read_densities %>% select(-width, -strand) %>%
    filter_at(4, all_vars(. > thr)) %>%
    na.omit()
  top_ranges = GRanges(
    seqnames = read_densities$seqnames,
    ranges = IRanges(start = read_densities$start,
                     end = read_densities$end)
  )
  return(top_ranges)
  
}

top_k27ac_over_cm =  function(k27ac_bigwig, top = 50) {
  name = strsplit(
    strsplit(
      k27ac_bigwig,
      "../data/bigwig/Yang_2019/GSM3314701_H3K27ac-"
    )[[1]][2],
    "_mm10.bw"
  )[[1]][1]
  
  read_densities = bw_loci(k27ac_bigwig,
                           loci = cm)
  read_densities = as_tibble(read_densities)
  thr = read_densities %>% pull(starts_with("GSM"))
  thr = thr[order(thr, decreasing = TRUE)][top]
  
  read_densities = read_densities %>% select(-width, -strand) %>%
    filter_at(4, all_vars(. > thr)) %>%
    na.omit()
  top_ranges = GRanges(
    seqnames = read_densities$seqnames,
    ranges = IRanges(start = read_densities$start,
                     end = read_densities$end)
  )
  return(top_ranges)
  
}

###### H3K27ac (Yang)
k27ac = yang[grep(pattern = "ac", x = yang)]
# execution
# loop through all pooled & scaled bigwig file
# lapply(epilc, function(x) {
#   lapply(k27ac, function(y) {
#     create_scatter(x, y)
#   })
# })

k27ac_correlations = lapply(epilc, function(x) {
  lapply(k27ac, function(y) {
    correlation_analysis(x, y)
  })
})

pearsons = bind_rows(k27ac_correlations)
pearsons = pearsons %>% filter(method == "pearson")
pearsons = pearsons %>% mutate(epilc_timepoint = case_when(str_detect(pearsons$comparison, "EpiLC_6h") ~ "6h",
                                                           str_detect(pearsons$comparison, "EpiLC_12h") ~ "12h",
                                                           str_detect(pearsons$comparison, "EpiLC_24h") ~ "24h",
                                                           str_detect(pearsons$comparison, "EpiLC_36h") ~ "36h",
                                                           str_detect(pearsons$comparison, "EpiLC_48h") ~ "48h",
                                                           str_detect(pearsons$comparison, "EpiLC_72h") ~ "72h"
                                                           ))
order = pearsons %>% arrange(desc(coefficient)) %>% pull(comparison)
order = factor(pearsons$comparison, levels = order)
fill = factor(pearsons$epilc_timepoint, levels = c("6h", "12h", "24h", "36h", "48h", "72h"))
ggplot(pearsons, aes(x = order, y = coefficient, fill = fill)) +
  geom_bar(stat="identity", color = "black") +
  ylim(0, 1) +
  scale_fill_manual(values = c("#fee0d2", "#fc9272", "#9ecae1", "#3182bd")) +
  theme_classic() +
  labs(
    title = "Correlations between H3.3 EpiLC and Yang K27ac ChIP-Seq over Cruz-Molina enhancers (+/- 1 kb)",
    x = "",
    y = "Pearson",
    fill = "EpiLC"
  ) +
  theme(
    text = element_text(size = 13),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 10, color = "black", angle = 90, vjust = 0.20, hjust = 1),
    axis.title.y = element_text(size = 15, color = "black"),
    axis.text.y = element_text(size = 15, color = "black")
  )

ggsave(
  glue("{result_folder}EpiLC_-K27ac_CruzMolina-correlations.png"),
  plot = last_plot(),
  width = 8,
  height = 5,
  dpi = 500
)
ggsave(
  glue("{result_folder}EpiLC_-K27ac_CruzMolina-correlations.pdf"),
  plot = last_plot(),
  width = 8,
  height = 5
)

pearsons = pearsons %>% mutate(epilc_timepoint_tipseq = case_when(str_detect(pearsons$comparison, "EpiLC_6h") ~ "EpiLC_6h",
                                                           str_detect(pearsons$comparison, "EpiLC_12h") ~ "EpiLC_12h",
                                                           str_detect(pearsons$comparison, "EpiLC_24h") ~ "EpiLC_24h",
                                                           str_detect(pearsons$comparison, "EpiLC_36h") ~ "EpiLC_36h",
                                                           str_detect(pearsons$comparison, "EpiLC_48h") ~ "EpiLC_48h",
                                                           str_detect(pearsons$comparison, "EpiLC_72h") ~ "EpiLC_72h"
))
pearsons = pearsons %>% mutate(epilc_timepoint_yang = case_when(str_detect(pearsons$comparison, "H3K27ac-0h") ~ "K27ac_0h",
                                                           str_detect(pearsons$comparison, "H3K27ac-1h") ~ "K27ac_1h",
                                                           str_detect(pearsons$comparison, "H3K27ac-6h") ~ "K27ac_6h",
                                                           str_detect(pearsons$comparison, "H3K27ac-12h") ~ "K27ac_12h",
                                                           str_detect(pearsons$comparison, "H3K27ac-24h") ~ "K27ac_24h",
                                                           str_detect(pearsons$comparison, "H3K27ac-36h") ~ "K27ac_36h",
                                                           str_detect(pearsons$comparison, "H3K27ac-48h") ~ "K27ac_48h",
                                                           str_detect(pearsons$comparison, "H3K27ac-72h") ~ "K27ac_72h"
))


# corr heatmap
mat = pivot_wider(pearsons, names_from = epilc_timepoint_tipseq, values_from = coefficient, id_cols = epilc_timepoint_yang) %>% 
  column_to_rownames(var = "epilc_timepoint_yang")
mat = as.matrix(mat)


png(
  file = glue("{result_folder}EpiLC_-K27ac_CruzMolina-correlations-pearsons_hm.png"),
  width = 13,
  height = 13,
  units = 'cm',
  res = 500
)
small_mat = mat[1:8, 1:4]
col_fun = colorRamp2(c(0, 0.5, 1), c("#9ecae1", "white", "#fc9272"))
Heatmap(
  mat,
  column_title = "H3.3 EpiLC vs. H3K27ac (Yang et al.) over mESC enhancers",
  row_title = "",
  name = "Pearson",
  clustering_method_rows = "complete",
  clustering_method_columns = "complete",
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 0.1),
  show_column_dend = TRUE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  heatmap_width = unit(10, "cm"),
  heatmap_height = unit(12, "cm"),
  row_names_gp = gpar(fontsize = 15),
  column_names_gp = gpar(fontsize = 12),
  column_title_gp = gpar(fontsize = 10),
  column_names_rot = 90,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", small_mat[i, j]), x, y, gp = gpar(fontsize = 10))
  }
)
dev.off()

###### H3K4me1 (Yang)
k4me1 = yang[grep(pattern = "H3K4me1", x = yang)]
# lapply(epilc, function(x) {
#   lapply(k4me1, function(y) {
#     create_scatter(x, y)
#   })
# })

k4me1_correlations = lapply(epilc, function(x) {
  lapply(k4me1, function(y) {
    correlation_analysis(x, y)
  })
})

pearsons = bind_rows(k4me1_correlations)
pearsons = pearsons %>% filter(method == "pearson")

pearsons = pearsons %>% mutate(epilc_timepoint_tipseq = case_when(str_detect(pearsons$comparison, "EpiLC_6h") ~ "EpiLC_6h",
                                                                  str_detect(pearsons$comparison, "EpiLC_12h") ~ "EpiLC_12h",
                                                                  str_detect(pearsons$comparison, "EpiLC_24h") ~ "EpiLC_24h",
                                                                  str_detect(pearsons$comparison, "EpiLC_36h") ~ "EpiLC_36h",
                                                                  str_detect(pearsons$comparison, "EpiLC_48h") ~ "EpiLC_48h",
                                                                  str_detect(pearsons$comparison, "EpiLC_72h") ~ "EpiLC_72h"
))
pearsons = pearsons %>% mutate(epilc_timepoint_yang = case_when(str_detect(pearsons$comparison, "H3K4me1-0h") ~ "K4me1_0h",
                                                                str_detect(pearsons$comparison, "H3K4me1-1h") ~ "K4me1_1h",
                                                                str_detect(pearsons$comparison, "H3K4me1-6h") ~ "K4me1_6h",
                                                                str_detect(pearsons$comparison, "H3K4me1-12h") ~ "K4me1_12h",
                                                                str_detect(pearsons$comparison, "H3K4me1-24h") ~ "K4me1_24h",
                                                                str_detect(pearsons$comparison, "H3K4me1-36h") ~ "K4me1_36h",
                                                                str_detect(pearsons$comparison, "H3K4me1-48h") ~ "K4me1_48h",
                                                                str_detect(pearsons$comparison, "H3K4me1-72h") ~ "K4me1_72h"
))


# corr heatmap
mat = pivot_wider(pearsons, names_from = epilc_timepoint_tipseq, values_from = coefficient, id_cols = epilc_timepoint_yang) %>% 
  column_to_rownames(var = "epilc_timepoint_yang")
mat = as.matrix(mat)

png(
  file = glue("{result_folder}EpiLC_-K4me1_CruzMolina-correlations-pearsons_hm.png"),
  width = 13,
  height = 13,
  units = 'cm',
  res = 500
)
small_mat = mat[1:8, 1:4]
col_fun = colorRamp2(c(0.5, 0.75, 1), c("#9ecae1", "white", "#fc9272"))
Heatmap(
  mat,
  column_title = "H3.3 EpiLC vs. H3K4me1 (Yang et al.) over mESC enhancers",
  row_title = "",
  name = "Pearson",
  clustering_method_rows = "complete",
  clustering_method_columns = "complete",
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 0.1),
  show_column_dend = TRUE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  heatmap_width = unit(10, "cm"),
  heatmap_height = unit(12, "cm"),
  row_names_gp = gpar(fontsize = 15),
  column_names_gp = gpar(fontsize = 12),
  column_title_gp = gpar(fontsize = 10),
  column_names_rot = 90,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", small_mat[i, j]), x, y, gp = gpar(fontsize = 10))
  }
)
dev.off()

# arrange scatter plots
pdfs = list.files(glue("{result_folder}"), pattern = "pdf", full.names = TRUE)
pdfs = pdfs[grep("H3K4me1", pdfs)]
sc_6h = pdfs[grep("EpiLC_6h", pdfs)]
sc_6h = sc_6h[c(1, 3, 7, 2, 4, 5, 6, 8)]

plots = lapply(sc_6h, image_read_pdf)
plots = lapply(plots, function(x) {image_ggplot(image = x[1])})
ggarrange(plotlist=plots)

ggsave(
  glue("{result_folder}EpiLC_6h-K4me1_CruzMolina.png"),
  plot = last_plot(),
  width = 6,
  height = 6,
  dpi = 500
)
ggsave(
  glue("{result_folder}EpiLC_6h-K4me1_CruzMolina.pdf"),
  plot = last_plot(),
  width = 8,
  height = 8
)

sc_12h = pdfs[grep("EpiLC_12h", pdfs)]
sc_12h = sc_12h[c(1, 3, 7, 2, 4, 5, 6, 8)]

plots = lapply(sc_12h, image_read_pdf)
plots = lapply(plots, function(x) {image_ggplot(image = x[1])})
ggarrange(plotlist=plots)

ggsave(
  glue("{result_folder}EpiLC_12h-K4me1_CruzMolina.png"),
  plot = last_plot(),
  width = 8,
  height = 8,
  dpi = 500
)
ggsave(
  glue("{result_folder}EpiLC_12h-K4me1_CruzMolina.pdf"),
  plot = last_plot(),
  width = 8,
  height = 8
)


sc_24h = pdfs[grep("EpiLC_24h", pdfs)]
sc_24h = sc_24h[c(1, 3, 7, 2, 4, 5, 6, 8)]

plots = lapply(sc_24h, image_read_pdf)
plots = lapply(plots, function(x) {image_ggplot(image = x[1])})
ggarrange(plotlist=plots)

ggsave(
  glue("{result_folder}EpiLC_24h-K4me1_CruzMolina.png"),
  plot = last_plot(),
  width = 8,
  height = 8,
  dpi = 500
)
ggsave(
  glue("{result_folder}EpiLC_24h-K4me1_CruzMolina.pdf"),
  plot = last_plot(),
  width = 8,
  height = 8
)

sc_48h = pdfs[grep("EpiLC_48h", pdfs)]
sc_48h = sc_48h[c(1, 3, 7, 2, 4, 5, 6, 8)]

plots = lapply(sc_48h, image_read_pdf)
plots = lapply(plots, function(x) {image_ggplot(image = x[1])})
ggarrange(plotlist=plots)

ggsave(
  glue("{result_folder}EpiLC_48h-K4me1_CruzMolina.png"),
  plot = last_plot(),
  width = 8,
  height = 8,
  dpi = 500
)
ggsave(
  glue("{result_folder}EpiLC_48h-K4me1_CruzMolina.pdf"),
  plot = last_plot(),
  width = 8,
  height = 8
)

## heatmaps
# Heatmaps about enhancers with the most variable H3.3 signal
k27ac_read_dens = bw_loci(k27ac, loci = cm)
k27ac_read_dens = as_tibble(k27ac_read_dens)
k27ac_read_dens = k27ac_read_dens %>% 
  na.omit(k27ac_read_dens) %>% 
  rowwise() %>% 
  mutate(., variance = sd(c_across(c(starts_with("GSM")))),
         mean = mean(c_across(c(starts_with("GSM"))))) %>%  
  dplyr::select(-width, -strand) %>% 
  mutate(range_id = paste(seqnames, start, end, sep = "_"))

k4me_read_dens = bw_loci(k4me1, loci = cm)
k4me_read_dens = as_tibble(k4me_read_dens)
k4me_read_dens = k4me_read_dens %>% 
  na.omit(k4me_read_dens) %>% 
  rowwise() %>% 
  mutate(., variance = sd(c_across(c(starts_with("GSM")))),
         mean = mean(c_across(c(starts_with("GSM"))))) %>%  
  dplyr::select(-width, -strand) %>% 
  mutate(range_id = paste(seqnames, start, end, sep = "_"))

h33_read_dens = bw_loci(epilc, loci = cm)
h33_read_dens = as_tibble(h33_read_dens)
h33_read_dens = h33_read_dens %>% 
  na.omit(h33_read_dens) %>% 
  rowwise() %>% 
  mutate(., variance = sd(c_across(c(starts_with("EpiLC")))),
         mean = mean(c_across(c(starts_with("EpiLC")))))
top_var = h33_read_dens %>% arrange(desc(variance)) %>% 
  ungroup() %>% 
  top_n(., n = 50, wt = variance) %>% 
  dplyr::select(-width, -strand) %>% 
  mutate(range_id = paste(seqnames, start, end, sep = "_")) %>% 
  select(range_id, EpiLC_6h_pseudobulk_RPGC, EpiLC_12h_pseudobulk_RPGC, EpiLC_24h_pseudobulk_RPGC,
         EpiLC_48h_pseudobulk_RPGC)
range_order = top_var %>% pull(range_id)

top_var_k27ac = k27ac_read_dens %>% filter(range_id %in% top_var$range_id) %>% 
  select(range_id, GSM3314703_H3K27ac.6h_mm10, GSM3314704_H3K27ac.12h_mm10, GSM3314705_H3K27ac.24h_mm10,
         GSM3314707_H3K27ac.48h_mm10) 
top_var_k4me = k4me_read_dens %>% filter(range_id %in% top_var$range_id) %>% 
  select(range_id, GSM3314719_H3K4me1.6h_mm10, GSM3314720_H3K4me1.12h_mm10, GSM3314721_H3K4me1.24h_mm10,
         GSM3314723_H3K4me1.48h_mm10) 

top_var = top_var %>% column_to_rownames(var = "range_id") %>% as.matrix
log_top_var = log2(top_var)
colnames(log_top_var) = c("H3.3, 6h", "H3.3, 12h", "H3.3, 24h", "H3.3, 48h")
top_var_k27ac = top_var_k27ac %>% column_to_rownames(var = "range_id")
top_var_k27ac = top_var_k27ac[range_order, ]
top_var_k27ac = as.matrix(top_var_k27ac)
log_top_var_k27ac = log2(top_var_k27ac)
colnames(log_top_var_k27ac) = c("K27ac, 6h", "K27ac, 12h", "K27ac, 24h", "K27ac, 48h")

top_var_k4me = top_var_k4me %>% column_to_rownames(var = "range_id")
top_var_k4me = top_var_k4me[range_order, ]
top_var_k4me = as.matrix(top_var_k4me)
log_top_var_k4me = log2(top_var_k4me)
colnames(log_top_var_k4me) = c("K4me1, 6h", "K4me1, 12h", "K4me1, 24h", "K4me1, 48h")

library("RColorBrewer")
display.brewer.pal(n = 3, name = 'Purples')
brewer.pal(n = 3, name = "Purples")

col_fun1 = colorRamp2(c(0, 1, 2), brewer.pal(n = 3, name = "Purples"))
top_var_hm = Heatmap(
  log_top_var,
  column_title = "",
  row_title = "most variable H3.3 signals at enhancers",
  name = "log2 read count",
  # clustering_method_rows = "complete",
  # clustering_method_columns = "complete",
  col = col_fun1,
  rect_gp = gpar(col = "black", lwd = 0.1),
  show_column_dend = FALSE,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_row_dend = FALSE,
  show_row_names = FALSE,
  heatmap_width = unit(4, "cm"),
  heatmap_height = unit(15, "cm"),
  row_names_gp = gpar(fontsize = 15),
  column_names_gp = gpar(fontsize = 12),
  column_title_gp = gpar(fontsize = 10),
  column_names_rot = 90)

display.brewer.pal(n = 3, name = 'Blues')

col_fun2 = colorRamp2(c(0, 1, 2), brewer.pal(n = 3, name = "Blues"))
top_var_k27ac_hm = Heatmap(
  log_top_var_k27ac,
  column_title = "",
  row_title = "most variable H3K27ac signals at enhancers",
  name = "log2 read count",
  # clustering_method_rows = "complete",
  # clustering_method_columns = "complete",
  col = col_fun2,
  rect_gp = gpar(col = "black", lwd = 0.1),
  show_column_dend = FALSE,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_row_dend = FALSE,
  show_row_names = FALSE,
  heatmap_width = unit(3, "cm"),
  heatmap_height = unit(15, "cm"),
  row_names_gp = gpar(fontsize = 15),
  column_names_gp = gpar(fontsize = 12),
  column_title_gp = gpar(fontsize = 10),
  column_names_rot = 90)

display.brewer.pal(n = 3, name = 'Reds')

col_fun3 = colorRamp2(c(0, 0.25, 0.5), brewer.pal(n = 3, name = "Reds"))
top_var_k4me1_hm = Heatmap(
  log_top_var_k4me,
  column_title = "",
  row_title = "most variable H3K4me1 signals at enhancers",
  name = "log2 read count",
  # clustering_method_rows = "complete",
  # clustering_method_columns = "complete",
  col = col_fun3,
  rect_gp = gpar(col = "black", lwd = 0.1),
  show_column_dend = FALSE,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_row_dend = FALSE,
  show_row_names = FALSE,
  heatmap_width = unit(3, "cm"),
  heatmap_height = unit(15, "cm"),
  row_names_gp = gpar(fontsize = 15),
  column_names_gp = gpar(fontsize = 12),
  column_title_gp = gpar(fontsize = 10),
  column_names_rot = 90)

png(
  file = glue("{result_folder}most_variable_H3.3_enhancers-hm.png"),
  width = 14,
  height = 16,
  units = 'cm',
  res = 500
)
top_var_hm + top_var_k27ac_hm + top_var_k4me1_hm
dev.off()

# Heatmaps of top H3.3 ranges at Cruz-Molina enhancers
heatmaps_of_tops = function(h33_bigwig) {
  name = strsplit(
    strsplit(
      h33_bigwig,
      "../results/Seurat/cluster_bigwigs/H3.3_EpiLC_timepoints/"
    )[[1]][2],
    "_pseudobulk_RPGC.bigwig"
  )[[1]][1]
  
  print(name)
  
  # K27ac over enhancers with highest H3.3 signal
  k27ac_read_dens = bw_loci(k27ac, loci = top_h33_over_cm(h33_bigwig, top = 52))
  k27ac_read_dens = as_tibble(k27ac_read_dens)
  k27ac_read_dens = k27ac_read_dens %>% 
    na.omit(k27ac_read_dens) %>% 
    rowwise() %>% 
    mutate(., variance = sd(c_across(c(starts_with("GSM")))),
           mean = mean(c_across(c(starts_with("GSM"))))) %>%  
    dplyr::select(-width, -strand) %>% 
    mutate(range_id = paste(seqnames, start, end, sep = "_"))
  
  # K4me1 over enhancers with highest H3.3 signal
  k4me_read_dens = bw_loci(k4me1, loci = top_h33_over_cm(h33_bigwig, top = 52))
  k4me_read_dens = as_tibble(k4me_read_dens)
  k4me_read_dens = k4me_read_dens %>% 
    na.omit(k4me_read_dens) %>% 
    rowwise() %>% 
    mutate(., variance = sd(c_across(c(starts_with("GSM")))),
           mean = mean(c_across(c(starts_with("GSM"))))) %>%  
    dplyr::select(-width, -strand) %>% 
    mutate(range_id = paste(seqnames, start, end, sep = "_"))
  
  # EpiLC H3.3 read densities over enhancers with highest H3.3 signal
  h33_read_dens = bw_loci(epilc, loci = top_h33_over_cm(epilc[grep(name, epilc)], top = 52))
  h33_read_dens = as_tibble(h33_read_dens)
  h33_read_dens = h33_read_dens %>% 
    na.omit(h33_read_dens) %>% 
    rowwise() %>% 
    mutate(., variance = sd(c_across(c(starts_with("EpiLC")))),
           mean = mean(c_across(c(starts_with("EpiLC")))))
  
  # top 50 variable ranges based on SD
  top_var = h33_read_dens %>% arrange(desc(variance)) %>% 
    ungroup() %>% 
    top_n(., n = 50, wt = variance) %>% 
    dplyr::select(-width, -strand) %>% 
    mutate(range_id = paste(seqnames, start, end, sep = "_")) %>% 
    select(range_id, EpiLC_6h_pseudobulk_RPGC, EpiLC_12h_pseudobulk_RPGC, EpiLC_24h_pseudobulk_RPGC,
           EpiLC_48h_pseudobulk_RPGC)
  range_order = top_var %>% pull(range_id)
  
  top_var_k27ac = k27ac_read_dens %>% filter(range_id %in% top_var$range_id) %>% 
    select(range_id, GSM3314703_H3K27ac.6h_mm10, GSM3314704_H3K27ac.12h_mm10, GSM3314705_H3K27ac.24h_mm10,
           GSM3314707_H3K27ac.48h_mm10) 
  top_var_k4me = k4me_read_dens %>% filter(range_id %in% top_var$range_id) %>% 
    select(range_id, GSM3314719_H3K4me1.6h_mm10, GSM3314720_H3K4me1.12h_mm10, GSM3314721_H3K4me1.24h_mm10,
           GSM3314723_H3K4me1.48h_mm10) 
  
  # EpiLC H3.3 heatmap input
  top_var = top_var %>% column_to_rownames(var = "range_id") %>% as.matrix
  log_top_var = log2(top_var + 2)
  colnames(log_top_var) = c("H3.3, 6h", "H3.3, 12h", "H3.3, 24h", "H3.3, 48h")
  
  # K27ac heatmap input
  top_var_k27ac = top_var_k27ac %>% column_to_rownames(var = "range_id")
  top_var_k27ac = top_var_k27ac[range_order, ]
  top_var_k27ac = as.matrix(top_var_k27ac)
  log_top_var_k27ac = log2(top_var_k27ac + 2)
  colnames(log_top_var_k27ac) = c("K27ac, 6h", "K27ac, 12h", "K27ac, 24h", "K27ac, 48h")
  
  # K4me1 heatmap input
  top_var_k4me = top_var_k4me %>% column_to_rownames(var = "range_id")
  top_var_k4me = top_var_k4me[range_order, ]
  top_var_k4me = as.matrix(top_var_k4me)
  log_top_var_k4me = log2(top_var_k4me + 2)
  colnames(log_top_var_k4me) = c("K4me1, 6h", "K4me1, 12h", "K4me1, 24h", "K4me1, 48h")
  
  # exclude rows with "NA"
  ranges_wo_nas_k27ac = rownames(log_top_var_k27ac)[which(rownames(log_top_var_k27ac) != "NA")]
  ranges_wo_nas_k4me1 = rownames(log_top_var_k4me)[which(rownames(log_top_var_k4me) != "NA")]
  wo_nas = intersect(ranges_wo_nas_k27ac, ranges_wo_nas_k4me1)
  log_top_var_k4me = log_top_var_k4me[wo_nas,]
  log_top_var_k27ac = log_top_var_k27ac[wo_nas,]
  log_top_var = log_top_var[wo_nas,]
  
  # test
  if(all(rownames(log_top_var) == rownames(log_top_var_k4me))) {
    "rownames are the same"
  }
  
  # pearson heatmap
  pearsons_k27ac = numeric()
  for (i in seq(1, dim(log_top_var)[1])) {
    coef = cor(log_top_var[i, ], log_top_var_k27ac[i, ], method = "pearson")
    pearsons_k27ac = c(pearsons_k27ac, round(coef, 2))
  }
  pearsons_k4me1 = numeric()
  for (i in seq(1, dim(log_top_var)[1])) {
    coef = cor(log_top_var[i, ], log_top_var_k4me[i, ], method = "pearson")
    pearsons_k4me1 = c(pearsons_k4me1, round(coef, 2))
  }
  pearsons = tibble(k27ac_pearsons = pearsons_k27ac, k4me1_pearsons = pearsons_k4me1)
  pearsons = as.matrix(pearsons)
  rownames(pearsons) = rownames(log_top_var)
  colnames(pearsons) = c("K27ac", "K4me1")
  
  # heatmaps
  col_fun1 = colorRamp2(c(0, 2, 4), brewer.pal(n = 3, name = "Purples"))
  top_hm = Heatmap(
    log_top_var,
    column_title = "",
    row_title = glue("Highest {name} H3.3 signals at enhancers"),
    name = "H3.3 log2 read count",
    clustering_method_rows = "complete",
    # clustering_method_columns = "complete",
    col = col_fun1,
    rect_gp = gpar(col = "black", lwd = 0.1),
    show_column_dend = FALSE,
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    show_row_dend = TRUE,
    show_row_names = FALSE,
    heatmap_width = unit(4, "cm"),
    heatmap_height = unit(15, "cm"),
    row_names_gp = gpar(fontsize = 15),
    column_names_gp = gpar(fontsize = 12),
    column_title_gp = gpar(fontsize = 10),
    column_names_rot = 90)
  top_hm
  
  max = round(max(log_top_var_k27ac))
  min = round(min(log_top_var_k27ac))
  middle = mean(c(max, min))
  col_fun2 = colorRamp2(c(min, middle, max), brewer.pal(n = 3, name = "Blues"))
  top_k27ac_hm = Heatmap(
    log_top_var_k27ac,
    column_title = "",
    row_title = glue("Highest {name} H3.3 signals at enhancers"),
    name = "H3K27ac log2 read count",
    # clustering_method_rows = "complete",
    # clustering_method_columns = "complete",
    col = col_fun2,
    rect_gp = gpar(col = "black", lwd = 0.1),
    show_column_dend = FALSE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_row_dend = FALSE,
    show_row_names = FALSE,
    heatmap_width = unit(2, "cm"),
    heatmap_height = unit(15, "cm"),
    row_names_gp = gpar(fontsize = 15),
    column_names_gp = gpar(fontsize = 12),
    column_title_gp = gpar(fontsize = 10),
    column_names_rot = 90)
  top_k27ac_hm
  
  max = round(max(log_top_var_k4me))
  min = round(min(log_top_var_k4me))
  middle = mean(c(max, min))
  col_fun3 = colorRamp2(c(min, middle, max), brewer.pal(n = 3, name = "Greens"))
  top_k4me1_hm = Heatmap(
    log_top_var_k4me,
    column_title = "",
    row_title = glue("Highest {name} H3.3 signals at enhancers"),
    name = "H3K4me1 log2 read count",
    # clustering_method_rows = "complete",
    # clustering_method_columns = "complete",
    col = col_fun3,
    rect_gp = gpar(col = "black", lwd = 0.1),
    show_column_dend = FALSE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_row_dend = FALSE,
    show_row_names = FALSE,
    heatmap_width = unit(2, "cm"),
    heatmap_height = unit(15, "cm"),
    row_names_gp = gpar(fontsize = 15),
    column_names_gp = gpar(fontsize = 12),
    column_title_gp = gpar(fontsize = 10),
    column_names_rot = 90)
  top_k4me1_hm
  
  pearson_fun = colorRamp2(c(-1, 0, 1), c("#3182bd", "white", "#fc9272"))
  pearsons_hm = Heatmap(
    pearsons,
    column_title = "",
    row_title = glue("Pearson - between {name} and ChIP-Seq signals"),
    name = "Pearson",
    # clustering_method_rows = "complete",
    # clustering_method_columns = "complete",
    col = pearson_fun,
    rect_gp = gpar(col = "black", lwd = 0.1),
    show_column_dend = FALSE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_row_dend = FALSE,
    show_row_names = FALSE,
    heatmap_width = unit(1, "cm"),
    heatmap_height = unit(15, "cm"),
    row_names_gp = gpar(fontsize = 15),
    column_names_gp = gpar(fontsize = 12),
    column_title_gp = gpar(fontsize = 10),
    column_names_rot = 90)
  pearsons_hm
  
  # exporting
  png(
    file = glue("{result_folder}highest_{name}_H3.3_enhancers-hm.png"),
    width = 15,
    height = 16,
    units = 'cm',
    res = 500
  )
  hms = top_hm + top_k27ac_hm + top_k4me1_hm + pearsons_hm
  print(hms)
  dev.off()
  
  pdf(
    file = glue("{result_folder}highest_{name}_H3.3_enhancers-hm.pdf"),
    width = 7,
    height = 7
  )
  hms = top_hm + top_k27ac_hm + top_k4me1_hm + pearsons_hm
  print(hms)
  dev.off()
  
}

# loop through all H3.3 sciTIP-Seq bigwigs
lapply(epilc, heatmaps_of_tops)

# Heatmaps of top H3K27ac ranges at Cruz-Molina enhancers 
heatmaps_of_top_k27ac = function(k27ac_bigwig) {
  
  name = strsplit(strsplit(
    strsplit(
      k27ac_bigwig,
      "../data/bigwig/Yang_2019/GSM"
    )[[1]][2],
    "_H3K27ac-"
  )[[1]][2],
  "_mm10.bw")[[1]][1]
  
  print(name)
  if(!name %in% c("6h", "12h", "24h", "48h")) {
    return(print("Time point does not overlap with H3.3 EpiLC time points"))
  } 
  
  # K27ac over enhancers with highest H3.3 signal
  k27ac_read_dens = bw_loci(k27ac, loci = top_k27ac_over_cm(k27ac_bigwig, top = 52))
  k27ac_read_dens = as_tibble(k27ac_read_dens)
  k27ac_read_dens = k27ac_read_dens %>% 
    na.omit(k27ac_read_dens) %>% 
    rowwise() %>% 
    mutate(., variance = sd(c_across(c(starts_with("GSM")))),
           mean = mean(c_across(c(starts_with("GSM"))))) %>%  
    dplyr::select(-width, -strand) %>% 
    mutate(range_id = paste(seqnames, start, end, sep = "_"))
  
  # K4me1 over enhancers with highest H3.3 signal
  k4me1 = yang[grep(pattern = "H3K4me1", x = yang)]
  k4me1 = k4me1[c(3,4,5,7)]
  k4me_read_dens = bw_loci(k4me1, loci = top_k27ac_over_cm(k27ac_bigwig, top = 52))
  k4me_read_dens = as_tibble(k4me_read_dens)
  k4me_read_dens = k4me_read_dens %>% 
    na.omit(k4me_read_dens) %>% 
    rowwise() %>% 
    mutate(., variance = sd(c_across(c(starts_with("GSM")))),
           mean = mean(c_across(c(starts_with("GSM"))))) %>%  
    dplyr::select(-width, -strand) %>% 
    mutate(range_id = paste(seqnames, start, end, sep = "_"))
  
  # EpiLC H3.3 read densities over enhancers with highest H3.3 signal
  h33_read_dens = bw_loci(epilc, loci = top_k27ac_over_cm(k27ac_bigwig, top = 52))
  h33_read_dens = as_tibble(h33_read_dens)
  h33_read_dens = h33_read_dens %>% 
    na.omit(h33_read_dens) %>% 
    rowwise() %>% 
    mutate(., variance = sd(c_across(c(starts_with("EpiLC")))),
           mean = mean(c_across(c(starts_with("EpiLC"))))) %>% 
    mutate(range_id = paste(seqnames, start, end, sep = "_"))
  
  # top 50 variable ranges based on SD
  top_var = k27ac_read_dens %>% arrange(desc(variance)) %>% 
    ungroup() %>% 
    top_n(., n = 50, wt = variance) %>% 
    mutate(range_id = paste(seqnames, start, end, sep = "_")) %>% 
    select(range_id, GSM3314703_H3K27ac.6h_mm10, GSM3314704_H3K27ac.12h_mm10, GSM3314705_H3K27ac.24h_mm10,
           GSM3314707_H3K27ac.48h_mm10)
  range_order = top_var %>% pull(range_id)
  
  top_var_h33 = h33_read_dens %>% filter(range_id %in% top_var$range_id) %>% 
    select(range_id, EpiLC_6h_pseudobulk_RPGC, EpiLC_12h_pseudobulk_RPGC, EpiLC_24h_pseudobulk_RPGC,
           EpiLC_48h_pseudobulk_RPGC)
  top_var_k4me = k4me_read_dens %>% filter(range_id %in% top_var$range_id) %>% 
    select(range_id, GSM3314719_H3K4me1.6h_mm10, GSM3314720_H3K4me1.12h_mm10, GSM3314721_H3K4me1.24h_mm10,
           GSM3314723_H3K4me1.48h_mm10) 
  
  # EpiLC H3.3 heatmap input
  top_var = top_var %>% column_to_rownames(var = "range_id") %>% as.matrix
  log_top_var = log2(top_var + 2)
  colnames(log_top_var) = c("K27ac, 6h", "K27ac, 12h", "K27ac, 24h", "K27ac, 48h")
  
  # K27ac heatmap input
  top_var_h33 = top_var_h33 %>% column_to_rownames(var = "range_id")
  top_var_h33 = top_var_h33[range_order, ]
  top_var_h33 = as.matrix(top_var_h33)
  log_top_var_h33 = log2(top_var_h33 + 2)
  colnames(log_top_var_h33) = c("H3.3, 6h", "H3.3, 12h", "H3.3, 24h", "H3.3, 48h")
  
  # K4me1 heatmap input
  top_var_k4me = top_var_k4me %>% column_to_rownames(var = "range_id")
  top_var_k4me = top_var_k4me[range_order, ]
  top_var_k4me = as.matrix(top_var_k4me)
  log_top_var_k4me = log2(top_var_k4me + 2)
  colnames(log_top_var_k4me) = c("K4me1, 6h", "K4me1, 12h", "K4me1, 24h", "K4me1, 48h")
  
  # exclude rows with "NA"
  ranges_wo_nas_h33 = rownames(log_top_var_h33)[which(rownames(log_top_var_h33) != "NA")]
  ranges_wo_nas_k4me1 = rownames(log_top_var_k4me)[which(rownames(log_top_var_k4me) != "NA")]
  wo_nas = intersect(ranges_wo_nas_h33, ranges_wo_nas_k4me1)
  log_top_var_k4me = log_top_var_k4me[wo_nas,]
  log_top_var_h33 = log_top_var_h33[wo_nas,]
  log_top_var = log_top_var[wo_nas,]
  
  # test
  if(all(rownames(log_top_var) == rownames(log_top_var_k4me))) {
    "rownames are the same"
  }
  
  # pearson heatmap
  pearsons_k27ac = numeric()
  for (i in seq(1, dim(log_top_var)[1])) {
    coef = cor(log_top_var[i, ], log_top_var_h33[i, ], method = "pearson")
    pearsons_k27ac = c(pearsons_k27ac, round(coef, 2))
  }
  pearsons_k4me1 = numeric()
  for (i in seq(1, dim(log_top_var_h33)[1])) {
    coef = cor(log_top_var_h33[i, ], log_top_var_k4me[i, ], method = "pearson")
    pearsons_k4me1 = c(pearsons_k4me1, round(coef, 2))
  }
  pearsons = tibble(k27ac_pearsons = pearsons_k27ac, k4me1_pearsons = pearsons_k4me1)
  pearsons = as.matrix(pearsons)
  rownames(pearsons) = rownames(log_top_var)
  colnames(pearsons) = c("K27ac", "K4me1")
  
  # heatmaps
  col_fun1 = colorRamp2(c(0, 2, 4), brewer.pal(n = 3, name = "Purples"))
  top_hm = Heatmap(
    log_top_var_h33,
    column_title = "",
    row_title = glue("Highest {name} H3K27ac (Yang et al) signals at enhancers"),
    name = "H3.3 log2 read count",
    # clustering_method_rows = "complete",
    # clustering_method_columns = "complete",
    col = col_fun1,
    rect_gp = gpar(col = "black", lwd = 0.1),
    show_column_dend = FALSE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_row_dend = TRUE,
    show_row_names = FALSE,
    heatmap_width = unit(2, "cm"),
    heatmap_height = unit(15, "cm"),
    row_names_gp = gpar(fontsize = 15),
    column_names_gp = gpar(fontsize = 12),
    column_title_gp = gpar(fontsize = 10),
    column_names_rot = 90)
  top_hm
  
  max = round(max(log_top_var))
  min = round(min(log_top_var))
  middle = mean(c(max, min))
  col_fun2 = colorRamp2(c(min, middle, max), brewer.pal(n = 3, name = "Blues"))
  top_k27ac_hm = Heatmap(
    log_top_var,
    column_title = "",
    row_title = glue("Highest {name} H3K27ac (Yang et al) signals at enhancers"),
    name = "H3K27ac log2 read count",
    clustering_method_rows = "complete",
    #clustering_method_columns = "complete",
    col = col_fun2,
    rect_gp = gpar(col = "black", lwd = 0.1),
    show_column_dend = FALSE,
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    show_row_dend = TRUE,
    show_row_names = FALSE,
    heatmap_width = unit(4, "cm"),
    heatmap_height = unit(15, "cm"),
    row_names_gp = gpar(fontsize = 15),
    column_names_gp = gpar(fontsize = 12),
    column_title_gp = gpar(fontsize = 10),
    column_names_rot = 90)
  top_k27ac_hm
  
  max = round(max(log_top_var_k4me))
  min = round(min(log_top_var_k4me))
  middle = mean(c(max, min))
  col_fun3 = colorRamp2(c(min, middle, max), brewer.pal(n = 3, name = "Greens"))
  top_k4me1_hm = Heatmap(
    log_top_var_k4me,
    column_title = "",
    row_title = glue("Highest {name} H3K27ac (Yang et al) signals at enhancers"),
    name = "H3K4me1 log2 read count",
    # clustering_method_rows = "complete",
    # clustering_method_columns = "complete",
    col = col_fun3,
    rect_gp = gpar(col = "black", lwd = 0.1),
    show_column_dend = FALSE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_row_dend = FALSE,
    show_row_names = FALSE,
    heatmap_width = unit(2, "cm"),
    heatmap_height = unit(15, "cm"),
    row_names_gp = gpar(fontsize = 15),
    column_names_gp = gpar(fontsize = 12),
    column_title_gp = gpar(fontsize = 10),
    column_names_rot = 90)
  top_k4me1_hm
  
  pearson_fun = colorRamp2(c(-1, 0, 1), c("#3182bd", "white", "#fc9272"))
  pearsons_hm = Heatmap(
    pearsons,
    column_title = "",
    row_title = glue("Pearson - between {name} and ChIP-Seq signals"),
    name = "Pearson",
    # clustering_method_rows = "complete",
    # clustering_method_columns = "complete",
    col = pearson_fun,
    rect_gp = gpar(col = "black", lwd = 0.1),
    show_column_dend = FALSE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_row_dend = FALSE,
    show_row_names = FALSE,
    heatmap_width = unit(1, "cm"),
    heatmap_height = unit(15, "cm"),
    row_names_gp = gpar(fontsize = 15),
    column_names_gp = gpar(fontsize = 12),
    column_title_gp = gpar(fontsize = 10),
    column_names_rot = 90)
  pearsons_hm
  
  # exporting
  png(
    file = glue("{result_folder}highest_Yang_EpiLC_K27ac_{name}_H3.3_enhancers-hm.png"),
    width = 15,
    height = 16,
    units = 'cm',
    res = 500
  )
  hms = top_k27ac_hm + top_hm + top_k4me1_hm + pearsons_hm
  print(hms)
  dev.off()
  
  pdf(
    file = glue("{result_folder}highest_Yang_EpiLC_K27ac_{name}_H3.3_enhancers-hm.pdf"),
    width = 7,
    height = 7
  )
  hms = top_k27ac_hm + top_hm + top_k4me1_hm + pearsons_hm
  print(hms)
  dev.off()
  
}

# loop through all H3.3 sciTIP-Seq bigwigs
lapply(k27ac, heatmaps_of_top_k27ac)

overlaps_with_cm = lapply(epilc_peaks, function(x) {
  x = fread(x)
  x = GRanges(
    seqnames = x$V1,
    ranges = IRanges(
      start = x$V2,
      end = x$V3,
      names = x$V6
    )
  )
  epilc_peak_cm_ol = findOverlaps(x, cm_orig_size, type = "any")
  ol_perc = round(length(unique(queryHits(epilc_peak_cm_ol))) / length(unique(x)) * 100, 2)
  return(ol_perc)
})
peak_lengths = lapply(epilc_peaks, function(x) {
  x = fread(x)
  x = as_tibble(x)
  mean = x %>% mutate(diff = V3 - V2) %>% pull(diff) %>% mean
  sd = x %>% mutate(diff = V3 - V2) %>% pull(diff) %>% sd
  output = tibble(mean = mean, sd = sd)
  return(output)
})
peak_lengths = rbindlist(peak_lengths)
peak_lengths = peak_lengths %>% mutate(peak = epilc_peaks)



