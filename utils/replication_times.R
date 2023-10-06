# packages
suppressPackageStartupMessages({
  library("Seurat")
  library("Signac")
  library("glue")
  library("tidyverse")
  library("Matrix")
  library("GenomicFeatures")
  library("data.table")
  library("ggpubr")
  library("RColorBrewer")
})

set.seed(42)

# result folder
result_folder = "../results/Seurat/"

# replication time data
early = fread("../data/bed/RepliTime_early_mm10.bed")
mid = fread("../data/bed/RepliTime_mid_mm10.bed")
late = fread("../data/bed/RepliTime_late_mm10.bed")

# sciTIP-Seq H3.2 - columns: cells above 1000 reads overall 
raw = fread(file = "../data/20230316_H3.2/count_tables/20230316_H3.2_read_counts-cells_above_1000reads_5k.tsv")
rows = raw$range
raw = raw %>% dplyr::select(-range)
raw_sparse = as(as.matrix(raw), "sparseMatrix")
rownames(raw_sparse) = rows
rm(raw_sparse)

# create GenomicRanges for overlap analysis
# sciTIP-Seq H3.2
h32 = tibble(range = rows)
h32 = h32 %>% 
  separate(., range, sep = "_", into = c("V1", "V2", "V3")) %>% 
  mutate(V2 = as.numeric(V2), V3 = as.numeric(V3))
h32$type = "H3.2"
h32 = GRanges(
  seqnames = h32$V1,
  ranges = IRanges(
    start = h32$V2,
    end = h32$V3,
    names = h32$type,
  )
)

# replication time analysis
early$type = "early"
early = GRanges(
  seqnames = early$V1,
  ranges = IRanges(
    start = early$V2,
    end = early$V3,
    names = early$type,
  )
)
mid$type = "mid"
mid = GRanges(
  seqnames = mid$V1,
  ranges = IRanges(
    start = mid$V2,
    end = mid$V3,
    names = mid$type,
  )
)
late$type = "late"
late = GRanges(
  seqnames = late$V1,
  ranges = IRanges(
    start = late$V2,
    end = late$V3,
    names = late$type,
  )
)

# find overlaps with ranges from replication time analysis
## early regions
early_ol = findOverlaps(h32, early, type = "within", ignore.strand = FALSE)
early_5kbs = as.tibble(h32[queryHits(early_ol)])

early_rows = as.tibble(early)
early_rows = early_rows %>% mutate(range = paste(seqnames, start, end, sep = "_")) %>% pull(range)
early_ol = as.data.frame(early_ol)
early_rows = early_rows[unique(early_ol$subjectHits)]

# early replicating regions in H3.2 dataset
early_h32 = as_tibble(h32[unique(early_ol$queryHits)])
write_tsv(early_h32, "../data/bed/early_H3.2_regions.bed", col_names = FALSE)
early_h32 = early_h32 %>% mutate(range = paste(seqnames, start, end, sep = "-")) %>% pull(range)
early_regions = raw[unique(early_ol$subjectHits),]

## mid regions
mid_ol = findOverlaps(h32, mid, type = "within", ignore.strand = TRUE)
mid_rows = as.tibble(mid)
mid_rows = mid_rows %>% mutate(range = paste(seqnames, start, end, sep = "_")) %>% pull(range)
mid_ol = as.data.frame(mid_ol)
mid_rows = mid_rows[unique(mid_ol$subjectHits)]

# mid replicating regions in H3.2 dataset
mid_h32 = as_tibble(h32[unique(mid_ol$queryHits)])
write_tsv(mid_h32, "../data/bed/mid_H3.2_regions.bed", col_names = FALSE)
mid_h32 = mid_h32 %>% mutate(range = paste(seqnames, start, end, sep = "-")) %>% pull(range)
mid_regions = raw[unique(mid_ol$subjectHits),]

## late regions
late_ol = findOverlaps(h32, late, type = "within", ignore.strand = TRUE)
late_rows = as.tibble(late)
late_rows = late_rows %>% mutate(range = paste(seqnames, start, end, sep = "_")) %>% pull(range)
late_ol = as.data.frame(late_ol)
late_rows = late_rows[unique(late_ol$subjectHits)]

# late replicating regions in H3.2 dataset
late_h32 = as_tibble(h32[unique(late_ol$queryHits)])
write_tsv(late_h32, "../data/bed/late_H3.2_regions.bed", col_names = FALSE)
late_h32 = late_h32 %>% mutate(range = paste(seqnames, start, end, sep = "-")) %>% pull(range)
late_regions = raw[unique(late_ol$subjectHits),]

early_sparse = as(as.matrix(early_regions), "sparseMatrix")
rownames(early_sparse) = early_rows
early_seurat = CreateSeuratObject(counts = early_sparse, project = "sciTIP_Seq")
early_meta = early_seurat@meta.data %>% mutate(replication = "early")
early_seurat@meta.data = early_meta

mid_sparse = as(as.matrix(mid_regions), "sparseMatrix")
rownames(mid_sparse) = mid_rows
mid_seurat = CreateSeuratObject(counts = mid_sparse, project = "sciTIP_Seq")
mid_meta = mid_seurat@meta.data %>% mutate(replication = "mid")
mid_seurat@meta.data = mid_meta

late_sparse = as(as.matrix(late_regions), "sparseMatrix")
rownames(late_sparse) = late_rows
late_seurat = CreateSeuratObject(counts = late_sparse, project = "sciTIP_Seq")
late_meta = late_seurat@meta.data %>% mutate(replication = "late")
late_seurat@meta.data = late_meta

coembed = merge(x = early_seurat, y = c(mid_seurat, late_seurat))
meta = coembed@meta.data

# normalization
# coembed = RunTFIDF(coembed)
# coembed = FindTopFeatures(coembed, min.cutoff = 'q0')
# coembed = RunSVD(coembed)
# 
# # Non-linear dimension reduction and clustering
# coembed = RunUMAP(object = coembed,
#                  reduction = 'lsi',
#                  dims = 2:30)
# coembed = FindNeighbors(object = coembed,
#                        reduction = 'lsi',
#                        dims = 2:30)
# coembed = FindClusters(object = coembed,
#                       verbose = FALSE,
#                       resolution = 0.5,
#                       algorithm = 3)
# 
# dim = DimPlot(object = coembed, group.by = "replication", label = FALSE, pt.size = 2, label.size = 7) + 
#   scale_color_brewer(palette = "Set3") +
#   xlim(-10, 10) + 
#   ylim(-10, 10) + 
#   ggtitle("H3.2 sciTIP-Seq") +
#   theme(
#     text = element_text(size = 25),
#     plot.title = element_text(size = 20),
#     axis.text.x = element_text(size = 25, color = "black"),
#     axis.text.y = element_text(size = 25, color = "black")
#   )
# dim

# quality plots
nCount_violin_clusters = VlnPlot(coembed, group.by = "replication", features = "nCount_RNA", pt.size = 0.1) +
  scale_fill_brewer(palette = "Set3") +
  ggtitle("") +
  xlab("replication time") + 
  ylab("read count") +
  ylim(0, 8000) +
  scale_y_continuous(breaks = seq(0, 8000, 2000)) +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black", angle = 0),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoLegend() +
  stat_compare_means(label.y = 6000, label.x = 1.5) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", label.y = 5500)

nCount_violin_clusters

ggsave(
  glue("{result_folder}replc_time-read_count_violins.pdf"),
  plot = nCount_violin_clusters,
  width = 7,
  height = 5,
  device = "pdf"
)

ggsave(
  glue("{result_folder}replc_time-read_count_violins.png"),
  plot = nCount_violin_clusters,
  width = 7,
  height = 5,
  dpi = 500
)


# ggsave(
#   glue("{result_folder}replc_time-UMAP.pdf"),
#   plot = dim,
#   width = 7,
#   height = 5,
#   device = "pdf"
# )
# 
# ggsave(
#   glue("{result_folder}replc_time-UMAP.png"),
#   plot = dim,
#   width = 7,
#   height = 5,
#   dpi = 500
# )

# H3.2 sciTIP-Seq Seurat (5 kbp)
seurat_5kb = readRDS("../data/20230316_H3.2/count_tables/20230316_H3.2_read_counts-cells_above_1000reads_5k.Rds")

top_cells = function(seurat_object, time, sparse_matrix) {
  seurat_object = FindVariableFeatures(seurat_object, selection.method = "disp", nfeatures = 500)
  top_variables = unlist(lapply(seurat_object@assays$RNA@var.features, function(x) {
    gsub("-", "_", x)
  }))
  sums = colSums(sparse_matrix[top_variables,])
  thr = quantile(sums, .75)
  tops = names(sums[sums > thr])
  
  meta = seurat_5kb@meta.data
  cell_ids = rownames(meta)
  meta = meta %>% mutate(cell_id = cell_ids) %>% 
    mutate(!!paste0(time, "_repl_status") := ifelse(cell_id %in% tops, time, ""))
  return(meta)
}

seurat_5kb@meta.data = top_cells(seurat_object = early_seurat, time = "early", sparse_matrix = early_sparse)

# visualizations
DimPlot(
  seurat_5kb,
  pt.size = 2,
  label.size = 7,
  group.by = 'early_repl_status',
  repel = FALSE,
  order = "early",
  raster = FALSE
) +
  xlim(-15, 15) +
  ylim(-15, 15) +
  scale_colour_manual(values = c("red", "grey")) +
  ggtitle("sciTIP-Seq early (RepliTime, early)") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoAxes() + NoLegend()
  
seurat_5kb@meta.data = top_cells(seurat_object = late_seurat, time = "late", sparse_matrix = late_sparse)

# visualizations
DimPlot(
  seurat_5kb,
  pt.size = 2,
  label.size = 7,
  group.by = 'late_repl_status',
  repel = FALSE,
  order = "late",
  raster = FALSE
) +
  xlim(-15, 15) +
  ylim(-15, 15) +
  scale_colour_manual(values = c("red", "grey")) +
  ggtitle("sciTIP-Seq early (RepliTime, late)") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoAxes() + NoLegend()

meta = seurat_5kb@meta.data

# marker regions: cluster 2 VS cluster 1 (Seurat)
#markers = FindMarkers(seurat_5kb, ident.1 = 2, ident.2 = 1, test.use = "LR", logfc.threshold = 0)
#markers = markers %>% mutate(range = rownames(markers))
#write_tsv(markers, glue("{result_folder}FindMarker-cluster2_vs_cluster1_LR.tsv"))
markers = read_tsv(glue("{result_folder}FindMarker-cluster2_vs_cluster1_LR.tsv"))

# overlaps of cluster 2 marker regions
cluster2_markers = markers %>% dplyr::filter(avg_log2FC > 0.2)

cluster2_markers = cluster2_markers %>% 
  separate(., range, sep = "-", into = c("V1", "V2", "V3")) %>% 
  mutate(V2 = as.numeric(V2), V3 = as.numeric(V3))
cluster2_markers$type = "cluster2_marker_region"
cluster2_markers = GRanges(
  seqnames = cluster2_markers$V1,
  ranges = IRanges(
    start = cluster2_markers$V2,
    end = cluster2_markers$V3,
    names = cluster2_markers$type,
  )
)
cl2_early = findOverlaps(cluster2_markers, early, type = "within", ignore.strand = FALSE)
length(cl2_early)
cl2_mid = findOverlaps(cluster2_markers, mid, type = "within", ignore.strand = FALSE)
length(cl2_mid)
cl2_late = findOverlaps(cluster2_markers, late, type = "within", ignore.strand = FALSE)
length(cl2_late)

cl2_overlaps = tibble(overlap = c(length(cl2_early), length(cl2_mid), length(cl2_late)), "time" = c("early", "mid", "late"), "Seurat_cluster" = "cluster 2")

# overlaps of cluster 1 marker regions 
cluster1_markers = markers %>% dplyr::filter(avg_log2FC < -0.2)

cluster1_markers = cluster1_markers %>% 
  separate(., range, sep = "-", into = c("V1", "V2", "V3")) %>% 
  mutate(V2 = as.numeric(V2), V3 = as.numeric(V3))
cluster1_markers$type = "cluster1_marker_region"
cluster1_markers = GRanges(
  seqnames = cluster1_markers$V1,
  ranges = IRanges(
    start = cluster1_markers$V2,
    end = cluster1_markers$V3,
    names = cluster1_markers$type,
  )
)
cl1_early = findOverlaps(cluster1_markers, early, type = "within", ignore.strand = FALSE)
length(cl1_early)
cl1_mid = findOverlaps(cluster1_markers, mid, type = "within", ignore.strand = FALSE)
length(cl1_mid)
cl1_late = findOverlaps(cluster1_markers, late, type = "within", ignore.strand = FALSE)
length(cl1_late)

cl1_overlaps = tibble(overlap = c(length(cl1_early), length(cl1_mid), length(cl1_late)), "time" = c("early", "mid", "late"), "Seurat_cluster" = "cluster 1")

# barplot about overlaps
overlaps = rbind(cl1_overlaps, cl2_overlaps)

order = factor(overlaps$time, levels = c("early", "mid", "late"))
ggplot(overlaps, aes(fill = order, y = overlap, x = Seurat_cluster)) + 
  geom_bar(position="dodge", stat="identity", color = "black") +
  scale_fill_manual(values = c("#bdbdbd", "#fec44f", "#9ecae1")) +
  labs(
    title = "Overlaps of marker regions (|log2FC| > 0.2) with RepliTime",
    x = "",
    y = "# of overlaps",
    fill = "S phase"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 13),
    plot.title = element_text(size = 11),
    axis.text.x = element_text(size = 23, color = "black", angle = 0),
    axis.title.y = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 23, color = "black")
  ) 
  
ggsave(
  glue("{result_folder}Marker_regions-replc_time-bars.pdf"),
  plot = last_plot(),
  width = 6,
  height = 4,
  device = "pdf"
)

ggsave(
  glue("{result_folder}Marker_regions-replc_time-bars.png"),
  plot = last_plot(),
  width = 6,
  height = 4,
  dpi = 500
)

top_cluster1 = markers %>% arrange(avg_log2FC) %>% top_n(-1, wt = avg_log2FC) %>% pull(range)
featureplot1 = FeaturePlot(seurat_5kb, features = top_cluster1) +
  xlim(-15, 15) +
  ylim(-15, 15) +
  scale_color_gradientn(colors = c("grey", "orange", "#c44a46")) +
  ggtitle("cluster 1 marker region") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 11),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoAxes() + NoLegend()
top_cluster2 = markers %>% arrange(avg_log2FC) %>% top_n(1, wt = avg_log2FC) %>% pull(range)
featureplot2 = FeaturePlot(seurat_5kb, features = top_cluster2) +
  xlim(-15, 15) +
  ylim(-15, 15) +
  scale_color_gradientn(colors = c("grey", "orange", "#c44a46")) +
  ggtitle("cluster 2 marker region") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 11),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoAxes() + NoLegend()

ggarrange(featureplot1, featureplot2)

ggsave(
  glue("{result_folder}Marker_regions-cl2_cl1.pdf"),
  plot = last_plot(),
  width = 4,
  height = 4,
  device = "pdf"
)

ggsave(
  glue("{result_folder}Marker_regions-cl2_cl1.png"),
  plot = last_plot(),
  width = 4,
  height = 4,
  dpi = 500
)


