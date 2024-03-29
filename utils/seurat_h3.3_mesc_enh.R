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
  library("ComplexHeatmap")
  library("circlize")
  library("ggpubr")
})

set.seed(42)

# result folder
result_folder = "../results/Seurat/"

# create Seurat object
raw = fread(file = "../data/20230510_EpiLC/count_tables/mESC_cm_only_enhancers-read_counts.tsv")
raw = raw %>% mutate(range = paste(raw$`#'chr'`, raw$`'start'`, raw$`'end'`, sep = "_"))
rows = raw$range
raw = raw %>% dplyr::select(-range, -`#'chr'`, -`'start'`, -`'end'`) 
cols = lapply(colnames(raw), function(x) {
  x = strsplit(x, split = "_L001_R1_001.fastq_mm10.q10.fixmate.psort.bam")[[1]][1]
  return(x)
})
cols = unlist(cols)
cols = lapply(cols, function(x) {
  x = gsub("'", "", x = x)
  return(x)
})
cols = unlist(cols)
colnames(raw) = cols
raw_sparse = as(as.matrix(raw), "sparseMatrix")
rownames(raw_sparse) = rows
rownames(raw) = rows

enh_ranked = rownames(raw)[order(rowSums(raw), decreasing = TRUE)]
head(enh_ranked)

seurat = CreateSeuratObject(counts = raw_sparse, project = "EpiLC_enhancers")

# export Rds
saveRDS(seurat, "../data/20230316_H3.2/count_tables/mESC_only_enhancers-H3.3.Rds")

# load Seurat object
seurat = readRDS(file = "../data/20230316_H3.2/count_tables/mESC_only_enhancers-H3.3.Rds")

# normalization
seurat = RunTFIDF(seurat)
seurat = FindTopFeatures(seurat, min.cutoff = 'q0')
seurat = RunSVD(seurat)

# Non-linear dimension reduction and clustering
seurat = RunUMAP(object = seurat,
                 reduction = 'lsi',
                 dims = 2:30)
seurat = FindNeighbors(object = seurat,
                       reduction = 'lsi',
                       dims = 2:30)
seurat = FindClusters(object = seurat,
                      verbose = FALSE,
                      resolution = 0.5,
                      algorithm = 3)

# quality plots
my_comparison = list(c("0", "1"))
nCount_violin_clusters = VlnPlot(seurat, group.by = "seurat_clusters", features = "nCount_RNA", pt.size = 0.1) +
  scale_fill_brewer(palette = "Set3") +
  ggtitle("Number of H3.3 counts over mESC enhancers") +
  scale_y_continuous(limits = c(0, 8000), breaks = seq(0, 8000, 2000)) +
  xlab("cluster") + 
  ylab("read count") +
  #ylim(0, 60000) +
  #scale_y_continuous(breaks = seq(0, 1000, 200)) +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black", angle = 0),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoLegend() +
  stat_compare_means(comparisons = my_comparison, label = "p.signif", method = "t.test",
                     ref.group = ".all.", label.y = 7400, tip.length = 0.05,
                     bracket.nudge.y = -12, step.increase = 0.8)

nCount_violin_clusters

ggsave(
  glue("{result_folder}nCount_over_mESC_enh-violins.png"),
  plot = nCount_violin_clusters,
  width = 7,
  height = 7,
  dpi = 300,
)

ggsave(
  glue("{result_folder}nCount_over_mESC_enh-violins.pdf"),
  plot = nCount_violin_clusters,
  width = 7,
  height = 7,
  device = "pdf"
)

cluster1 = seurat@meta.data %>% filter(seurat_clusters == 1)
cluster1 = tibble("cells_with_H3.3_over_mESC_enhancers" = rownames(cluster1))
write_tsv(cluster1, glue("{result_folder}cells_with_H3.3_over_mESC_enhancers.tsv"))
cluster1_highcount = seurat@meta.data %>% filter(seurat_clusters == 1) %>% 
  filter(nCount_RNA >= quantile(nCount_RNA, .75))
cluster1_highcount = tibble("cells_with_highH3.3_over_mESC_enhancers" = rownames(cluster1_highcount))
write_tsv(cluster1_highcount, glue("{result_folder}cells_with_highH3.3_over_mESC_enhancers.tsv"))

# visualizations
dim = DimPlot(object = seurat, label = FALSE, pt.size = 2, label.size = 7) + 
  scale_color_brewer(palette = "Set3") +
  xlim(-20, 20) + 
  ylim(-20, 20) + 
  ggtitle("H3.3 sciTIP-Seq") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )
dim

# add sample information
meta = seurat@meta.data
meta = meta %>% mutate(sample_id = as.numeric(unname(sapply(rownames(meta), function(x) strsplit(x, split = "_")[[1]][3]))))
meta = meta %>% 
  mutate(sample = case_when((sample_id >= 1 & sample_id <= 96) ~ "EpiLC_6h",
                            (sample_id > 96 & sample_id <= 192) ~ "EpiLC_12h",
                            (sample_id > 192 & sample_id <= 288) ~ "EpiLC_24h",
                            (sample_id > 288 & sample_id <= 384) ~ "EpiLC_48h", 
                            TRUE ~ "-")) %>% dplyr::select(-sample_id)
seurat@meta.data = meta

seurat@meta.data$sample <- factor(seurat@meta.data$sample, levels = c("EpiLC_6h", "EpiLC_12h", "EpiLC_24h", "EpiLC_48h"))

# quality plots
nCount_violin_samples = VlnPlot(seurat, group.by = "sample", features = "nCount_RNA", pt.size = 0.1) +
  scale_fill_brewer(palette = "Set3") +
  ggtitle("read counts at mESC enhancers") +
  xlab("sample") + 
  ylab("read count") + 
  #ylim(0, 100) +
  scale_y_continuous(limits = c(0, 8000), breaks = seq(0, 8000, 2000)) +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 14, color = "black", angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoLegend() +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", label.y = 7500, tip.length = 0.05,
                     bracket.nudge.y = -12, step.increase = 0.8)
nCount_violin_samples

dim_samples = DimPlot(object = seurat, label = FALSE, pt.size = 2, group.by = "sample") + 
  scale_color_brewer(palette = "Set3") +
  xlim(-10, 10) + 
  ylim(-10, 10) + 
  ggtitle("") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )
dim_samples

# ggsave(
#   glue("{result_folder}Seurat_H3.3_sciTIP_0510-sample_level_UMAP.png"),
#   plot = dim_samples,
#   width = 10,
#   height = 10,
#   dpi = 300,
# )
# 
# ggsave(
#   glue("{result_folder}Seurat_H3.3_sciTIP_0510-sample_level_UMAP.pdf"),
#   plot = dim_samples,
#   width = 10,
#   height = 10,
#   device = "pdf"
# )

dims = ggarrange(dim, dim_samples)

ggsave(
  glue("{result_folder}Seurat_H3.3_sciTIP_0510-UMAPs.png"),
  plot = dims,
  width = 12,
  height = 5,
  dpi = 300,
)

ggsave(
  glue("{result_folder}Seurat_H3.3_sciTIP_0510-UMAPs.pdf"),
  plot = dims,
  width = 12,
  height = 5,
  device = "pdf"
)

violins = ggarrange(nCount_violin_clusters, nCount_violin_samples)

ggsave(
  glue("{result_folder}Seurat_H3.3_sciTIP_0510-quality_plots-mESCenhancers.png"),
  plot = violins,
  width = 12,
  height = 5,
  dpi = 300,
)

ggsave(
  glue("{result_folder}Seurat_H3.3_sciTIP_0510-quality_plots-mESCenhancers.pdf"),
  plot = violins,
  width = 12,
  height = 5,
  device = "pdf"
)

# H3.3 marker regions
unique(seurat@meta.data$sample)

# at level of time points
marker_analysis = FindMarkers(seurat, test.use = "roc",
                              ident.1 = "EpiLC_6h",
                              ident.2 = "EpiLC_48h",
                              group.by = "sample",
                              logfc.threshold = 0,
                              min.pct = 0)
# summary = marker_analysis %>% 
#   dplyr::filter(abs(avg_log2FC) > 0.1) %>% 
#   dplyr::filter(p_val_adj < 0.05)

# on whole dataset
all_marker_analysis = FindAllMarkers(
  seurat,
  test.use = "roc",
  group.by = "sample",
  logfc.threshold = 0,
  min.pct = 0
)
