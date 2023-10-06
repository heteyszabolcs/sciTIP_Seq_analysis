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
               "circlize") 


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
                           loci = ol_w_cm,)
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

