if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse",
               "data.table",
               "wigglescout",
               "ggpubr",
               "glue",
               "GenomicRanges") 

# output folder
result_folder = "../results/ChIP-Seq/"

peaks = "../data/ChIP-Seq/Buecker_EpiLC_GSE56138/MACS2/"
peaks = list.files(peaks, pattern = ".narrowPeak", full.names = TRUE)

# Cruz-Molina active enhancers (mESC)
cm = "../data/bed/ESC_Enhancer_CruzMolina.active_mm10.bed"
cm = fread(cm)
cm = cm %>% mutate(V2 = V2 - 500, V3 = V3 + 500)

cm$V6 = "Cruz-Molina_active_enh"
cm = GRanges(
  seqnames = cm$V1,
  ranges = IRanges(
    start = cm$V2,
    end = cm$V3,
    names = cm$V6
  )
)

strong_peaks = function(peak_set, name) {
  peak_set = fread(peak_set)
  above_med = peak_set %>% filter(V7 > median(peak_set$V7))
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

oct4 = peaks[grep("p3", peaks)]
oct4_rep1 = strong_peaks(oct4[1], name = "EpiLC_oct4_rep1")
oct4_rep2 = strong_peaks(oct4[2], name = "EpiLC_oct4_rep2")
ol = findOverlaps(oct4_rep1, oct4_rep2, type = "any")
oct4 = oct4_rep1[queryHits(ol)]

enh = Reduce(subsetByOverlaps, list(k27ac, p300, oct4))
length(enh)
head(enh)
enh_int = Reduce(intersect, list(k27ac, p300, oct4))
length(enh_int)
head(enh_int)

epilc_enh_bed = as.data.frame(enh_int)
epilc_enh_bed = epilc_enh_bed[,c(1,2,3)]

write_tsv(epilc_enh_bed, glue("{result_folder}GSE56138_EpiLC_enhancers-Oct4_p300_K27ac.bed"), 
          col_names = FALSE)

ol = findOverlaps(enh_int, cm, type = "any")
epilc_only_enh = enh_int[-queryHits(ol)]
length(epilc_only_enh)
head(epilc_only_enh)

epilc_enh_bed = as.data.frame(epilc_only_enh)
epilc_enh_bed = epilc_enh_bed[,c(1,2,3)]
write_tsv(epilc_enh_bed, glue("{result_folder}GSE56138_EpiLC_only_enhancers-Oct4_p300_K27ac.bed"), 
          col_names = FALSE)

cm_only_enh = cm[-subjectHits(ol)]

cm_enh_bed = as_tibble(cm_only_enh)
cm_enh_bed = cm_enh_bed[,c(1,2,3)]
write_tsv(cm_enh_bed, glue("{result_folder}ESC_Enhancer_CruzMolina_only_enhancers.bed"), 
          col_names = FALSE)
