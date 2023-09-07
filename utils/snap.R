# packages
suppressPackageStartupMessages({
  library("Seurat")
  library("Signac")
  library("glue")
  library("tidyverse")
  library("Matrix")
  library("GenomicFeatures")
  library("GenomicRanges")
  library("data.table")
  library("ggpubr")
  library("purrr")
  library("SnapATAC")
  library("umap")
})

set.seed(42)

# result folder
result_folder = "../results/Seurat/"
# Seurat cluster ids
cluster_ids = list.files("../results/Seurat/", pattern = "*_ids", full.names = TRUE)

# load Seurat object + marker analysis output
seurat =
  readRDS(file = "../data/20230316_H3.2/count_tables/20230316_H3.2_read_counts-cells_above_1000reads.Rds")
markers = fread("../results/Seurat/Seurat_H3.2_sciTIP_Marker_analysis.tsv")
markers = markers %>% 
  dplyr::filter(avg_log2FC > 0.2) %>% 
  dplyr::filter(p_val_adj < 0.05)

mat = as.matrix(seurat@assays$sciTIP_Seq_H3.2@counts)
bins = tibble("V1" = rownames(mat)) %>% separate("V1", into = c("V1", "V2", "V3")) %>%
  mutate(V4 = "test",
         V2 = as.numeric(V2),
         V3 = as.numeric(V3))
bins = GRanges(
  seqnames = bins$V1,
  ranges = IRanges(
    start = bins$V2,
    end = bins$V3,
    names = bins$V4,
  )
)

## generate snap file for H3.2 sciTip-Seq
snap = createSnapFromBmat(as(t(as.matrix(mat)), "sparseMatrix"),
                          barcodes = colnames(mat),
                          bins = bins)

chr.exclude = seqlevels(snap@feature)[grep("random|chrM", seqlevels(snap@feature))]

idy = grep(paste(chr.exclude, collapse = "|"), snap@feature)

if (length(idy) > 0) {
  snap = snap[, -idy, mat = "bmat"]
}

bin.cov = log10(Matrix::colSums(snap@bmat) + 1)
bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95)
idy = which(bin.cov <= bin.cutoff & bin.cov > 0)
snap = snap[, idy, mat = "bmat"]

# keep bins with read count above 1000
idx = which(Matrix::rowSums(snap@bmat) > 1000)
snap = snap[idx, ]

# binarize sciTIP-Seq data
snap = makeBinary(snap, mat = "bmat")

# scRepli-seq integration
screpliseq_mids = list.files("../data/GSE108556_scRepli-Seq/midS_phase/", patter = "binary.bedGraph", 
                             full.names = TRUE)
screpliseq_earlys = list.files("../data/GSE108556_scRepli-Seq/early_25percentile_S_phase/", patter = "binary.bedGraph", 
                             full.names = TRUE)
screpliseq_lates = list.files("../data/GSE108556_scRepli-Seq/late_75percentile_S_phase/", patter = "binary.bedGraph", 
                               full.names = TRUE)

# return consistently replicating (scRepli-seq signal = 1) regions
retrieve_repl_regions = function(screpliseq_datasets) {
  datasets = lapply(screpliseq_datasets, fread)
  datasets = datasets %>% purrr::reduce(inner_join, by = c("V1", "V2", "V3"))
  cols = sapply(seq(screpliseq_datasets), function(x)
    paste0("cell_", as.character(x)))
  colnames(datasets) = c("chr", "start", "end", cols)
  replicating_regions = datasets %>% dplyr::filter(if_all(-c(chr, start, end), ~ . == 1))
  replicating_regions$type = "midS"
  replicating_regions = GRanges(
    seqnames = replicating_regions$chr,
    ranges = IRanges(
      start = replicating_regions$start,
      end = replicating_regions$end,
      names = replicating_regions$type,
    )
  )
  return(replicating_regions)
}

# retrieve replicating sites for all S subphase
screpliseq_mids = retrieve_repl_regions(screpliseq_datasets = screpliseq_mids)
screpliseq_earlys = retrieve_repl_regions(screpliseq_datasets = screpliseq_earlys)
screpliseq_lates = retrieve_repl_regions(screpliseq_datasets = screpliseq_lates)

# keep cluster specific cells
retrieve_clusters = function(cluster_ids, snap) {
  cluster = fread(cluster_ids, header = FALSE)
  cluster = sapply(cluster$V1, function(x) {
    paste0(x, ".fastq")
  })
  cluster = intersect(rownames(snap@bmat), cluster)
  cluster = snap@bmat[cluster, ]
  cluster = tibble(region = colnames(cluster)) %>%
    separate(region, sep = "-", into = c("V1", "V2", "V3")) %>% mutate(V2 = as.numeric(V2),
                                                                       V3 = as.numeric(V3))
  cluster$type = "h3.2_bins"
  cluster = GRanges(
    seqnames = cluster$V1,
    ranges = IRanges(
      start = cluster$V2,
      end = cluster$V3,
      names = cluster$type,
    )
  )
  
  return(cluster)
  
}

# scRepli-seq (80 kb bins) overlap with sciTip-seq (5 kb bins)
overlaps = function(cluster, screpliseq) {
  ol = findOverlaps(cluster, screpliseq, type = "any", ignore.strand = FALSE)
  both = cluster[queryHits(ol)]
  ids_both = tibble(seq = as.vector(both@seqnames), start = both@ranges@start,
                    end = as.numeric(start) + 5000, merged = paste(seq, start, end, sep = "-")) %>% 
    distinct(merged) %>% pull(merged)  
  return(ids_both)
}

mod_cluster_ids = function(cluster_ids) {
  cluster = fread(cluster_ids, header = FALSE)
  cluster = sapply(cluster$V1, function(x) {
    paste0(x, ".fastq")
  })
  return(cluster)
}

# Seurat cluster 0 of sciTIP-Seq over regions replicating in mid-S phase
screpliseq_sets = list(screpliseq_earlys, screpliseq_mids, screpliseq_lates)

subsetting = function(screpliseq_set, cluster, cluster_ids) {
  subset = snap@bmat[intersect(mod_cluster_ids(cluster_ids), rownames(snap@bmat)),
                     intersect(overlaps(cluster = cluster, screpliseq = screpliseq_set), 
                               colnames(snap@bmat))]
  
}

# input for visualizations
creating_hist_input = function(subset, cluster_label) {
  early = tibble(total_rc = rowSums(subset[[1]]), cluster = cluster_label, S_phase = "early")
  mid = tibble(total_rc = rowSums(subset[[2]]), cluster = cluster_label, S_phase = "mid")
  late = tibble(total_rc = rowSums(subset[[3]]), cluster = cluster_label, S_phase = "late")
  hist_input = rbind(early, mid, late)
}

# H3.2 visualization of replicating sites in early/mid/late S phase cells
# histogram
hist = function(hist_input, title) {
  ggplot(hist_input, aes(x = total_rc, fill = S_phase)) +
    geom_histogram(position = "identity", alpha = 0.4) +
    geom_density(alpha = 1.0) +
    scale_fill_manual(values = c("#fc9272", "#a6bddb", "#99d8c9")) +
    ylim(0,400) +
    labs(title = title,
         x = "",
         y = "Density",
         fill = "S phase") +
    theme_classic() +
    theme(
      text = element_text(size = 20),
      plot.title = element_text(size = 15),
      axis.text.x = element_text(size = 8, color = "black"))
}

# boxplot
# y axis: sum of H3.2 binary signal from the snap file
bp = function(hist_input, title) {
  hist_input = hist_input %>% dplyr::filter(total_rc < quantile(total_rc, 0.95))
  order = factor(hist_input$S_phase, levels = c("early", "mid", "late"))
  ggplot(hist_input, aes(y = log10(total_rc), x = order, fill = S_phase)) +
    geom_boxplot(color = "black") +
    scale_fill_manual(values = c("#fc9272", "#a6bddb", "#99d8c9")) +
    ylim(0, 10) +
    labs(title = title,
         x = "",
         y = "log10(total signal)",
         fill = "S phase") +
    theme_classic() +
    theme(legend.position = "none") +
    theme(
      text = element_text(size = 20),
      plot.title = element_text(size = 15),
      axis.text.x = element_text(size = 8, color = "black")) +
    stat_compare_means(label.y = 9,
                       label.x = 1.25,
                       size = 2)
}
  
# create subsets
# cluster specific H3.2 regions VS. early/late/mid S phase replicating sites
subsets_cluster0 = lapply(screpliseq_sets, subsetting,
                          cluster = retrieve_clusters(cluster_ids = cluster_ids[1], snap = snap), 
                          cluster_ids = cluster_ids[1])
subsets_cluster1 = lapply(screpliseq_sets, subsetting, 
                          cluster = retrieve_clusters(cluster_ids = cluster_ids[2], snap = snap),
                          cluster_ids = cluster_ids[2])
subsets_cluster2 = lapply(screpliseq_sets, subsetting,
                          cluster = retrieve_clusters(cluster_ids = cluster_ids[3], snap = snap),
                          cluster_ids = cluster_ids[3])
subsets_cluster3 = lapply(screpliseq_sets, subsetting, 
                          cluster = retrieve_clusters(cluster_ids = cluster_ids[4], snap = snap),
                          cluster_ids = cluster_ids[4])
subsets_cluster4 = lapply(screpliseq_sets, subsetting,
                          cluster = retrieve_clusters(cluster_ids = cluster_ids[5], snap = snap),
                          cluster_ids = cluster_ids[5])
subsets_cluster5 = lapply(screpliseq_sets, subsetting, 
                          cluster = retrieve_clusters(cluster_ids = cluster_ids[6], snap = snap),
                          cluster_ids = cluster_ids[6])

input_cluster0 = creating_hist_input(subset = subsets_cluster0, cluster_label = "0")
input_cluster1 = creating_hist_input(subset = subsets_cluster1, cluster_label = "1")
input_cluster2 = creating_hist_input(subset = subsets_cluster2, cluster_label = "2")
input_cluster3 = creating_hist_input(subset = subsets_cluster3, cluster_label = "3")
input_cluster4 = creating_hist_input(subset = subsets_cluster4, cluster_label = "4")
input_cluster5 = creating_hist_input(subset = subsets_cluster5, cluster_label = "5")

hist0 = hist(hist_input = input_cluster0, title = "cluster 0")
hist1 = hist(hist_input = input_cluster1, title = "cluster 1")
hist2 = hist(hist_input = input_cluster2, title = "cluster 2")
hist3 = hist(hist_input = input_cluster3, title = "cluster 3")
hist4 = hist(hist_input = input_cluster4, title = "cluster 4")
hist5 = hist(hist_input = input_cluster5, title = "cluster 5")
(hist0 + hist1 + hist2) / (hist3 + hist4 + hist5)

ggsave(
  glue("{result_folder}H3.2_sciTIP-Seq-scRepli-Seq_cluster_read_counts_hist.png"),
  plot = last_plot(),
  width = 12,
  height = 6,
  dpi = 300,
)

ggsave(
  glue("{result_folder}H3.2_sciTIP-Seq-scRepli-Seq_cluster_read_counts_hist.pdf"),
  plot = last_plot(),
  width = 12,
  height = 6
)

bp0 = bp(hist_input = input_cluster0, title = "cluster 0")
bp1 = bp(hist_input = input_cluster1, title = "cluster 1")
bp2 = bp(hist_input = input_cluster2, title = "cluster 2")
bp3 = bp(hist_input = input_cluster3, title = "cluster 3")
bp4 = bp(hist_input = input_cluster4, title = "cluster 4")
bp5 = bp(hist_input = input_cluster5, title = "cluster 5")
(bp0 + bp1 + bp2) / (bp3 + bp4 + bp5)

ggsave(
  glue("{result_folder}H3.2_sciTIP-Seq-scRepli-Seq_cluster_read_counts_box.png"),
  plot = last_plot(),
  width = 12,
  height = 6,
  dpi = 300,
)

ggsave(
  glue("{result_folder}H3.2_sciTIP-Seq-scRepli-Seq_cluster_read_counts_box.pdf"),
  plot = last_plot(),
  width = 12,
  height = 6
)

compute_marker_regions_gr = function(cluster_label) {
  cluster_markers = markers %>% dplyr::filter(cluster == as.numeric(cluster_label)) %>% 
    separate(gene, into = c("chr", "start", "end"), sep = "-")  %>% 
    mutate(start = as.numeric(start), end = as.numeric(end))
  cluster_markers$type = glue("Seurat_cluster{as.character(cluster_label)}")
  cluster_markers = GRanges(
    seqnames = cluster_markers$chr,
    ranges = IRanges(
      start = cluster_markers$start,
      end = cluster_markers$end,
      names = cluster_markers$type,
    )
  )
}

compute_marker_regions = function(cluster) {
  cluster_markers = markers %>% dplyr::filter(cluster == as.numeric(cluster))
  return(cluster_markers)
}

screpliseq_mids = list.files("../data/GSE108556_scRepli-Seq/midS_phase/", 
                             patter = "binary.bedGraph", 
                             full.names = TRUE)
screpliseq_earlys = list.files("../data/GSE108556_scRepli-Seq/early_25percentile_S_phase/", 
                               patter = "binary.bedGraph", 
                               full.names = TRUE)
screpliseq_lates = list.files("../data/GSE108556_scRepli-Seq/late_75percentile_S_phase/", 
                              patter = "binary.bedGraph", 
                              full.names = TRUE)

replication_state = function(screpliseq_sets,
                             marker_regions) {
  screpliseq_set = fread(screpliseq_sets[1])
  screpliseq_set$type = "scRepli-seq"
  screpliseq_set_gr = GRanges(
    seqnames = screpliseq_set$V1,
    ranges = IRanges(
      start = screpliseq_set$V2,
      end = screpliseq_set$V3,
      names = screpliseq_set$type,
    )
  )
  ol = findOverlaps(screpliseq_set_gr,
                    marker_regions,
                    type = "any",
                    ignore.strand = FALSE)
  both = screpliseq_set_gr[queryHits(ol)]
  
  repl_signals = lapply(screpliseq_sets, function(x) {
    df = fread(x)
    repl_signal = as_tibble(both) %>% distinct_all() %>% dplyr::select(-width, -strand) %>%
      inner_join(., df, by = c(
        "seqnames" = "V1",
        "start" = "V2",
        "end" = "V3"
      )) %>% pull("V4")
  })
  
  repl_matrix = matrix(nrow = length(repl_signals[[1]]),
                       ncol = length(repl_signals))
  for (i in seq(length(repl_signals))) {
    repl_matrix[, i] = repl_signals[[i]]
  }
  return(repl_matrix)
}

make_heatmap = function(replication_state_binary, scirepliseq_label, cluster_label) {
  library("ComplexHeatmap")
  library("circlize")
  col_fun = colorRamp2(c(-1, 0, 1), c("black", "white", "#fec44f"))
  lgd3 = Legend(labels = c("repl.", "no data", "non-repl."), 
                legend_gp = gpar(fill = 7:9, color = c("black", "white", "#fec44f")), 
                title = "state")
  hm = Heatmap(
    replication_state_binary,
    column_title = glue("{scirepliseq_label} scRepli-seq"),
    row_title = glue("H3.2 sciTIP-Seq, {cluster_label}"),
    name = "replication state",
    row_km = 2,
    column_km = 2,
    #clustering_method_rows = "complete",
    col = col_fun,
    #rect_gp = gpar(col = "black", lwd = 0.1),
    #top_annotation = ha,
    show_column_dend = FALSE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_row_dend = FALSE,
    heatmap_width = unit(3, "cm"),
    heatmap_height = unit(6, "cm"),
    row_names_gp = gpar(fontsize = 2),
    column_names_gp = gpar(fontsize = 6),
    column_names_rot = 90,
    show_heatmap_legend = FALSE
  )
  draw(hm)
  draw(lgd3, x = unit(0.70, "npc"), y = unit(0.5, "npc"))
}


## visualize the replication states of the clusters based on scRepli-Seq data
# sciTIP-Seq cluster 1
cluster1_early = replication_state(screpliseq_sets = screpliseq_earlys,
                             marker_regions = compute_marker_regions_gr(cluster = 1))
cluster1_mid = replication_state(screpliseq_sets = screpliseq_mids,
                             marker_regions = compute_marker_regions_gr(cluster = 1))
cluster1_late = replication_state(screpliseq_sets = screpliseq_lates,
                                 marker_regions = compute_marker_regions_gr(cluster = 1))

pdf(
  file = glue("{result_folder}cluster1_early_sciRepli-Seq.pdf"),
  width = 6,
  height = 6
)
cl1_early_hm = make_heatmap(replication_state_binary = cluster1_early,
             scirepliseq_label = "early", 
             cluster_label = "cluster 1")
dev.off()
pdf(
  file = glue("{result_folder}cluster1_mid_sciRepli-Seq.pdf"),
  width = 6,
  height = 6
)
cl1_mid_hm = make_heatmap(replication_state_binary = cluster1_mid,
             scirepliseq_label = "mid", 
             cluster_label = "cluster 1")
dev.off()
pdf(
  file = glue("{result_folder}cluster1_late_sciRepli-Seq.pdf"),
  width = 6,
  height = 6
)
cl1_late_hm = make_heatmap(replication_state_binary = cluster1_late,
             scirepliseq_label = "late", 
             cluster_label = "cluster 1")
dev.off()
cols = c(
  "scRNA-Seq" = "#bdbdbd",
  "0" = "#bdbdbd",
  "1" = "#de2d26",
  "2" = "#bdbdbd",
  "3" = "#bdbdbd",
  "4" = "#bdbdbd",
  "5" = "#bdbdbd")

DimPlot(
  seurat,
  pt.size = 2,
  label.size = 7,
  group.by = 'seurat_clusters',
  repel = TRUE,
  order = "3",
  raster = TRUE
) +
  xlim(-15, 15) +
  ylim(-15, 15) +
  scale_colour_manual(values = cols) +
  ggtitle("sciTIP-Seq cluster 1") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoAxes() + NoLegend()

ggsave(
  glue("{result_folder}Seurat_H3.2_sciTIP_UMAP_cluster1.pdf"),
  plot = last_plot(),
  width = 4,
  height = 4
)

# sciTIP-Seq cluster 2
cluster2_early = replication_state(screpliseq_sets = screpliseq_earlys,
                                   marker_regions = compute_marker_regions_gr(cluster = 2))
cluster2_mid = replication_state(screpliseq_sets = screpliseq_mids,
                                 marker_regions = compute_marker_regions_gr(cluster = 2))
cluster2_late = replication_state(screpliseq_sets = screpliseq_lates,
                                  marker_regions = compute_marker_regions_gr(cluster = 2))

pdf(
  file = glue("{result_folder}cluster2_early_sciRepli-Seq.pdf"),
  width = 6,
  height = 6
)
cl2_early_hm = make_heatmap(replication_state_binary = cluster2_early,
                            scirepliseq_label = "early", 
                            cluster_label = "cluster 2")
dev.off()
pdf(
  file = glue("{result_folder}cluster2_mid_sciRepli-Seq.pdf"),
  width = 6,
  height = 6
)
cl2_mid_hm = make_heatmap(replication_state_binary = cluster2_mid,
                          scirepliseq_label = "mid", 
                          cluster_label = "cluster 2")
dev.off()
pdf(
  file = glue("{result_folder}cluster2_late_sciRepli-Seq.pdf"),
  width = 6,
  height = 6
)
cl2_late_hm = make_heatmap(replication_state_binary = cluster2_late,
                           scirepliseq_label = "late", 
                           cluster_label = "cluster 2")
dev.off()
cols = c(
  "scRNA-Seq" = "#bdbdbd",
  "0" = "#bdbdbd",
  "1" = "#bdbdbd",
  "2" = "#de2d26",
  "3" = "#bdbdbd",
  "4" = "#bdbdbd",
  "5" = "#bdbdbd")

DimPlot(
  seurat,
  pt.size = 2,
  label.size = 7,
  group.by = 'seurat_clusters',
  repel = TRUE,
  order = "3",
  raster = TRUE
) +
  xlim(-15, 15) +
  ylim(-15, 15) +
  scale_colour_manual(values = cols) +
  ggtitle("sciTIP-Seq cluster 2") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoAxes() + NoLegend()

ggsave(
  glue("{result_folder}Seurat_H3.2_sciTIP_UMAP_cluster2.pdf"),
  plot = last_plot(),
  width = 4,
  height = 4
)

# sciTIP-Seq cluster 3
cluster3_early = replication_state(screpliseq_sets = screpliseq_earlys,
                                   marker_regions = compute_marker_regions_gr(cluster = 3))
cluster3_mid = replication_state(screpliseq_sets = screpliseq_mids,
                                 marker_regions = compute_marker_regions_gr(cluster = 3))
cluster3_late = replication_state(screpliseq_sets = screpliseq_lates,
                                  marker_regions = compute_marker_regions_gr(cluster = 3))

pdf(
  file = glue("{result_folder}cluster3_early_sciRepli-Seq.pdf"),
  width = 6,
  height = 6
)
cl3_early_hm = make_heatmap(replication_state_binary = cluster3_early,
                            scirepliseq_label = "early", 
                            cluster_label = "cluster 3")
dev.off()
pdf(
  file = glue("{result_folder}cluster3_mid_sciRepli-Seq.pdf"),
  width = 6,
  height = 6
)
cl3_mid_hm = make_heatmap(replication_state_binary = cluster3_mid,
                          scirepliseq_label = "mid", 
                          cluster_label = "cluster 3")
dev.off()
pdf(
  file = glue("{result_folder}cluster3_late_sciRepli-Seq.pdf"),
  width = 6,
  height = 6
)
cl3_late_hm = make_heatmap(replication_state_binary = cluster3_late,
                           scirepliseq_label = "late", 
                           cluster_label = "cluster 3")
dev.off()
cols = c(
  "scRNA-Seq" = "#bdbdbd",
  "0" = "#bdbdbd",
  "1" = "#bdbdbd",
  "2" = "#bdbdbd",
  "3" = "#de2d26",
  "4" = "#bdbdbd",
  "5" = "#bdbdbd")

DimPlot(
  seurat,
  pt.size = 2,
  label.size = 7,
  group.by = 'seurat_clusters',
  repel = TRUE,
  order = "3",
  raster = TRUE
) +
  xlim(-15, 15) +
  ylim(-15, 15) +
  scale_colour_manual(values = cols) +
  ggtitle("sciTIP-Seq cluster 3") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoAxes() + NoLegend()

ggsave(
  glue("{result_folder}Seurat_H3.2_sciTIP_UMAP_cluster3.pdf"),
  plot = last_plot(),
  width = 4,
  height = 4
)