library("data.table")
library("GenomicRanges")
library("tidyverse")
library("wigglescout")

peaks_6h = fread("../results/SEACR/EpiLC_6h.stringent_thr0.2.bed")
peaks_12h = fread("../results/SEACR/EpiLC_12h.stringent_thr0.2.bed")
peaks_24h = fread("../results/SEACR/EpiLC_24h.stringent_thr0.2.bed")
peaks_48h = fread("../results/SEACR/EpiLC_48h.stringent_thr0.2.bed")

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

decreasing_H33 = read_counts[which(apply(mat, MARGIN = 1, monotonic_decreas)),]



    