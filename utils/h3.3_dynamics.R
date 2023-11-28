# packages
suppressPackageStartupMessages({
  library("glue")
  library("tidyverse")
  library("GenomicFeatures")
  library("data.table")
  library("ggpubr")
  library("wigglescout")
  library("ComplexHeatmap")
  library("circlize")
})

# output folder
result_folder = "../results/wigglescout/"

# EpiLC H3.3 sciTIP-Seq peaks
epilc_peaks = list.files("../results/SEACR/",
                         pattern = "*bed",
                         full.names = TRUE)

# EpiLC H3.3 sciTIP-Seq signals
epilc = list.files("../results/Seurat/cluster_bigwigs/H3.3_EpiLC_timepoints/",
                   pattern = "*bigwig",
                   full.names = TRUE)

# mESC H3.3 signal
h33 = "../data/ChIP-Seq/mESC_GSE59189/mESC_H3.3_L001_R1_001_RPGC.bigwig"

# get top H3.3 peaks
top_peaks = function(peak) {
  bed = fread(peak)
  thr = quantile(bed$V5, .75)
  bed = bed %>% dplyr::filter(V5 > thr)
  return(bed)
}
epilc_peaks_top = lapply(epilc_peaks, top_peaks)
epilc_peaks_top = rbindlist(epilc_peaks_top)
epilc_peaks_top = epilc_peaks_top %>% group_by(V6) %>% summarize(max(V5)) %>% 
  separate(V6, sep = ":", into = c("chr", "range")) %>% 
  separate(range, sep = "-", into = c("start", "end"), remove = FALSE) %>% 
  mutate(start = as.numeric(start), end = as.numeric(end))

epilc_peaks_top$V6 = "EpiLC H3.3"
epilc_peaks_top = GRanges(
  seqnames = epilc_peaks_top$chr,
  ranges = IRanges(
    start = epilc_peaks_top$start,
    end = epilc_peaks_top$end,
    names = epilc_peaks_top$V6
  )
)

# 48 hour EpiLC enhancers: p300, Oct4, K27ac overlap
epilc_enh = fread("../data/bed/GSE56138_EpiLC_only_enhancers-Oct4_p300_K27ac.bed")
epilc_enh = epilc_enh %>% mutate(V2 = V2 - 1000, V3 = V3 + 1000)
epilc_enh$V6 = "EpiLC only enh"
epilc_enh = GRanges(
  seqnames = epilc_enh$V1,
  ranges = IRanges(
    start = epilc_enh$V2,
    end = epilc_enh$V3,
    names = epilc_enh$V6
  )
)

# overlap between EpiLC enhancers and highest EpiLC H3.3 peaks
ol = findOverlaps(epilc_peaks_top, epilc_enh)
epilc_h33_on_enh = epilc_peaks_top[queryHits(ol)]

# retrieve read densities by wigglescout and compute fold changes (to 6th hour)
read_densities = bw_loci(c(h33, epilc), loci = epilc_h33_on_enh)
read_densities = as.data.frame(read_densities)

reference = read_densities$mESC_H3.3_L001_R1_001_RPGC

fc_6h = log2((read_densities$EpiLC_6h_pseudobulk_RPGC + 1) / (reference + 1))
fc_12h = log2((read_densities$EpiLC_12h_pseudobulk_RPGC + 1) / (reference + 1))
fc_24h = log2((read_densities$EpiLC_24h_pseudobulk_RPGC + 1) / (reference + 1))
fc_48h = log2((read_densities$EpiLC_48h_pseudobulk_RPGC + 1) / (reference + 1))
fc_mat = tibble(seqnames = read_densities$seqnames, 
                start = read_densities$start,
                end = read_densities$end,
                "mESC" = read_densities$mESC_H3.3_L001_R1_001_RPGC,
                "epilc_6h" = fc_6h, "epilc_12h" = fc_12h, "epilc_24h" = fc_24h, "epilc_48h" = fc_48h)

# heatmap
hm_input = fc_mat[,4:ncol(fc_mat)]

pdf(
  file = glue("{result_folder}EpiLC-H3.3_dynamics_over_EpiLC_enhancers.pdf"),
  width = 4,
  height = 8
)
col_fun = colorRamp2(c(-4, 0, 4), c("#9ecae1", "white", "#fc9272"))
hm = Heatmap(
  as.matrix(hm_input[,2:5]),
  column_title = "",
  row_title = "H3.3 signal changes at EpiLC enhancers",
  name = "log2 fc",
  km = 5,
  col = col_fun,
  show_column_dend = TRUE,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  show_column_names = TRUE,
  show_row_names = FALSE,
  heatmap_width = unit(3, "cm"),
  heatmap_height = unit(15, "cm"))
hm
dev.off()

# k-means clustering
set.seed(42)
hm_input_scaled = scale(hm_input)
km = kmeans(hm_input_scaled, 5, nstart = 25)

clusters = list()
for (i in seq(1, 5)) {
  cluster = read_densities[which(km$cluster == i), ]
  cluster = cluster %>%
    dplyr::select(
      mESC_H3.3_L001_R1_001_RPGC,
      EpiLC_6h_pseudobulk_RPGC,
      EpiLC_12h_pseudobulk_RPGC,
      EpiLC_24h_pseudobulk_RPGC,
      EpiLC_48h_pseudobulk_RPGC
    ) %>%
    mutate(id = seq(1, nrow(.)))
  clusters[[i]] = cluster
}
clusters

# line plots
for(i in seq(1, 5)) {
  pdf(
    file = glue("{result_folder}EpiLC-H3.3_dyn_over_EpiLC_enh-cluster{i}.pdf"),
    width = 4,
    height = 4
  )
  print(matplot(
    t(clusters[[i]][, 1:4]),
    type = "b",
    pch = 1,
    ylab = glue("cluster {i}"),
    col = "#9ecae1",
    xaxt = 'n'
  ))
  means = colMeans(clusters[[i]][, 1:4])
  lines(means, col = alpha("#fc9272", 0.4), lwd = 12)
  axis(side = 1,at=c(1,2,3,4,5),labels=c("mESC", "6h","12h","24h","48h"))
  dev.off()
}

for(i in seq(1, 5)) {
  # cluster specific tibbles
  cluster = read_densities[which(km$cluster == i), ]
  cluster = cluster %>% dplyr::select(seqnames, start, end)
  write_tsv(cluster, glue("{result_folder}EpiLC-H3.3_dyn_over_EpiLC_enh-cluster{i}.tsv"),
            col_names = FALSE)
}

# cluster specific tibbles
cluster5 = read_densities[which(km$cluster == 5), ]
cluster5 = cluster5 %>% dplyr::select(EpiLC_6h_pseudobulk_RPGC, EpiLC_12h_pseudobulk_RPGC,
                                    EpiLC_24h_pseudobulk_RPGC, EpiLC_48h_pseudobulk_RPGC, everything())

cluster1 = read_densities[which(km$cluster == 1), ]
cluster1 = cluster5 %>% dplyr::select(EpiLC_6h_pseudobulk_RPGC, EpiLC_12h_pseudobulk_RPGC,
                                      EpiLC_24h_pseudobulk_RPGC, EpiLC_48h_pseudobulk_RPGC, everything())

# ordering bed files
# cluster 4
k27ac = "../data/ChIP-Seq/Buecker_EpiLC_GSE56138/H3K27ac_L001_R1_001_RPGC.bigwig"

cluster5 = fread("../results/wigglescout/EpiLC-H3.3_dyn_over_EpiLC_enh-cluster5.tsv")
cluster5$V4 = "EpiLC H3.3"
cluster5 = GRanges(
  seqnames = cluster5$V1,
  ranges = IRanges(
    start = cluster5$V2,
    end = cluster5$V3,
    names = cluster5$V4
  )
)

order = bw_loci(bwfiles = k27ac, loci = cluster5)
order = as.data.frame(order) 
order = order %>% arrange(desc(H3K27ac_L001_R1_001_RPGC))
order = order %>% dplyr::select(seqnames, start, end)
write_tsv(order, "../data/bed/kmeans_cluster5_ordered.bed", 
          col_names = FALSE)

# mESC H3.3 ChIP-Seq peaks - GSE59189
h33_peaks = fread("../data/ChIP-Seq/mESC_GSE59189/MACS2/mESC_H3.3_L001_R1_001_peaks.narrowPeak")
h33_peaks = h33_peaks %>% arrange(desc(V7)) %>% dplyr::select(V1, V2, V3)
h33_peaks = h33_peaks[1:100,]
write_tsv(h33_peaks, "../data/ChIP-Seq/mESC_GSE59189/MACS2/mESC_H3.3-top_100.bed", col_names = FALSE)



