if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(
  "tidyverse",
  "data.table",
  "wigglescout",
  "ggpubr",
  "glue",
  "GenomicRanges",
  "cotools"
)

# output folder
result_folder = "../results/ChIP-Seq/"

peaks = "../data/ChIP-Seq/Buecker_et_al_EpiLC_GSE56138/MACS2/"
peaks = list.files(peaks, pattern = ".narrowPeak", full.names = TRUE)

# Cruz-Molina active enhancers (mESC)
cm = "../data/bed/ESC_Enhancer_CruzMolina.active_mm10.bed"
cm = fread(cm)
length = cm %>% mutate(diff = V3 - V2) %>% pull(diff) %>% mean
sd = cm %>% mutate(diff = V3 - V2) %>% pull(diff) %>% sd
print(
  glue(
    "Mean length of Cruz-Molina enhancers is {as.character(round(length, 1))} 
    with SD {as.character(round(sd, 1))}"
  )
)
cm = cm %>% mutate(V2 = V2 - 500, V3 = V3 + 500)

cm$V6 = "Cruz-Molina_active_enh"
cm = GRanges(seqnames = cm$V1,
             ranges = IRanges(
               start = cm$V2,
               end = cm$V3,
               names = cm$V6
             ))

strong_peaks = function(peak_set, name) {
  peak_set = fread(peak_set)
  above_med = peak_set %>% filter(V7 > quantile(peak_set$V7, .25))
  above_med$V11 = name
  above_med = GRanges(
    seqnames = above_med$V1,
    ranges = IRanges(
      start = above_med$V2,
      end = above_med$V3,
      names = above_med$V11
    )
  )
}

# Buecker et al. EpiLC ChIP-Seq
k27ac = peaks[grep("K27ac", peaks)]
k27ac_rep1 = strong_peaks(k27ac[1], name = "EpiLC_K27ac_rep1")
k27ac_rep2 = strong_peaks(k27ac[2], name = "EpiLC_K27ac_rep2")
ol = findOverlaps(k27ac_rep1, k27ac_rep2, type = "any")
k27ac = k27ac_rep1[queryHits(ol)]

p300 = peaks[grep("p3", peaks)]
p300_rep1 = strong_peaks(p300[1], name = "EpiLC_p300_rep1")
p300_rep2 = strong_peaks(p300[2], name = "EpiLC_p300_rep2")
ol = findOverlaps(p300_rep1, p300_rep2, type = "any")
p300 = p300_rep1[queryHits(ol)]

oct4 = peaks[grep("Oct", peaks)]
oct4_rep1 = strong_peaks(oct4[1], name = "EpiLC_oct4_rep1")
oct4_rep2 = strong_peaks(oct4[2], name = "EpiLC_oct4_rep2")
ol = findOverlaps(oct4_rep1, oct4_rep2, type = "any")
oct4 = oct4_rep1[queryHits(ol)]

# overlaps
enh = Reduce(subsetByOverlaps, list(k27ac, p300, oct4))
length(enh)
head(enh)

# colocs of k27ac, p300, oct4
enh_int = Reduce(GenomicRanges::intersect, list(k27ac, p300, oct4))
length(enh_int)
head(enh_int)

epilc_enh_bed = as.data.frame(enh_int)
epilc_enh_bed = epilc_enh_bed[, c(1, 2, 3)]

write_tsv(
  epilc_enh_bed,
  glue(
    "{result_folder}GSE56138_EpiLC_enhancers-Oct4_p300_K27ac.bed"
  ),
  col_names = FALSE
)

# colocs of p300 and K27ac
enh_int2 = Reduce(GenomicRanges::intersect, list(k27ac, p300))
length(enh_int2)
head(enh_int2)

epilc_enh_bed2 = as.data.frame(enh_int2)
epilc_enh_bed2 = epilc_enh_bed2[, c(1, 2, 3)]

write_tsv(
  epilc_enh_bed2,
  glue(
    "{result_folder}GSE56138_EpiLC_enhancers-p300_K27ac.bed"
  ),
  col_names = FALSE
)

ol = findOverlaps(enh_int, cm, type = "any")
epilc_only_enh = enh_int[-queryHits(ol)]
length(epilc_only_enh)
head(epilc_only_enh)

epilc_enh_bed = as.data.frame(epilc_only_enh)
epilc_enh_bed = epilc_enh_bed[, c(1, 2, 3)]
length = epilc_enh_bed %>% mutate(diff = end - start) %>% pull(diff) %>% mean
sd = epilc_enh_bed %>% mutate(diff = end - start) %>% pull(diff) %>% sd
print(
  glue(
    "Mean length of EpiLC enhancers is {as.character(round(length, 1))} 
    with SD {as.character(round(sd, 1))}"
  )
)
write_tsv(
  epilc_enh_bed,
  glue(
    "{result_folder}GSE56138_EpiLC_only_enhancers-Oct4_p300_K27ac.bed"
  ),
  col_names = FALSE
)

ol = findOverlaps(enh_int2, cm, type = "any")
epilc_only_enh2 = enh_int2[-queryHits(ol)]
length(epilc_only_enh2)
head(epilc_only_enh2)

epilc_enh_bed = as.data.frame(epilc_only_enh2)
epilc_enh_bed = epilc_enh_bed[, c(1, 2, 3)]

write_tsv(
  epilc_enh_bed,
  glue(
    "{result_folder}GSE56138_EpiLC_only_enhancers-p300_K27ac.bed"
  ),
  col_names = FALSE
)

cm_only_enh = cm[-subjectHits(ol)]

cm_enh_bed = as_tibble(cm_only_enh)
cm_enh_bed = cm_enh_bed[, c(1, 2, 3)]
write_tsv(
  cm_enh_bed,
  glue(
    "{result_folder}ESC_Enhancer_CruzMolina_only_enhancers.bed"
  ),
  col_names = FALSE
)

# quantify overlaps
ol = findOverlaps(enh_int, cm, type = "any")
both = length(queryHits(ol))
cm_only = length(cm[-subjectHits(ol)])
epilc_only = length(enh_int[-queryHits(ol)])

### Jaccard analysis with stringent H3.3 peaks
# strong h3.3 sciTIP-Seq peaks
h33_peaks_path = "../results/SEACR/"
h33_peaks = list.files(h33_peaks_path, pattern = "*.bed")

create_gr_object = function(bed, name) {
  bed = fread(glue("{h33_peaks_path}{bed}"))
  bed$type = name
  bed = GRanges(
    seqnames = bed$V1,
    ranges = IRanges(
      start = bed$V2,
      end = bed$V3,
      names = bed$type
    )
  )
}

h33_peaks_6h = create_gr_object(h33_peaks[grep("6h", x = h33_peaks)], name = "6h")
h33_peaks_12h = create_gr_object(h33_peaks[grep("12h", x = h33_peaks)], name = "12h")
h33_peaks_24h = create_gr_object(h33_peaks[grep("24h", x = h33_peaks)], name = "24h")
h33_peaks_48h = create_gr_object(h33_peaks[grep("48h", x = h33_peaks)], name = "48h")

peaks = list(
  h33_peaks_6h,
  h33_peaks_12h,
  h33_peaks_24h,
  h33_peaks_48h,
  k27ac,
  oct4,
  p300,
  cm,
  epilc_only_enh
)
jaccards = matrix(NA_real_, length(peaks), length(peaks))
colnames(jaccards) = c(
  "EpiLC 6h",
  "EpiLC 12h",
  "EpiLC 24h",
  "EpiLC 48h",
  "EpiLC H3K27ac",
  "EpiLC Oct4",
  "EpiLC p300",
  "mESC enhancers",
  "EpiLC enhancers"
)
rownames(jaccards) = c(
  "EpiLC 6h",
  "EpiLC 12h",
  "EpiLC 24h",
  "EpiLC 48h",
  "EpiLC H3K27ac",
  "EpiLC Oct4",
  "EpiLC p300",
  "mESC enhancers",
  "EpiLC enhancers"
)
for (i in seq(1, ncol(jaccards))) {
  for (j in seq(1, nrow(jaccards))) {
    jaccard = genomicCorr.jaccard(peaks[[i]], peaks[[j]])
    jaccards[i, j] = jaccard
  }
}
is.matrix(jaccards)

pdf(file = "../results/genomic_ranges/EpiLC_H3.3_SEACR_peaks-jaccard_matrix.pdf",
    width = 7,
    height = 6)
col_fun = colorRamp2(c(0, 0.25, 0.5), c("#9ecae1", "white", "#fc9272"))
Heatmap(
  jaccards,
  column_title = "",
  row_title = "",
  name = "Jaccard index",
  # row_km = 2,
  # column_km = 1,
  clustering_method_rows = "complete",
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 0.1),
  #top_annotation = ha,
  show_column_dend = TRUE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  heatmap_width = unit(14, "cm"),
  heatmap_height = unit(14, "cm"),
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 90,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", jaccards[i, j]), x, y, gp = gpar(fontsize = 10))
  }
)
dev.off()
