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
raw = fread(file = "../data/20230510_EpiLC/count_tables/EpiLC_only_enhancers-k27ac_cluster1-read_counts.tsv")
raw = raw %>% mutate(range = paste(raw$`#'chr'`, raw$`'start'`, raw$`'end'`, sep = "_"))
rows = raw$range
raw = raw %>% dplyr::select(-range, -`#'chr'`, -`'start'`, -`'end'`) 
cols = lapply(colnames(raw), function(x) {
  x = strsplit(x, split = "_L001_R1_001.fastq_mm10.q10.fixmate.psort.bam")[[1]][1]
  return(x)
})
cols = unlist(cols)
colnames(raw) = cols
raw = as.matrix(raw)
rownames(raw) = rows
raw = raw[,colSums(raw) > 10]
raw = raw[rowSums(raw) > 0,]
raw_sparse = as(as.matrix(raw), "sparseMatrix")
rownames(raw_sparse) = rownames(raw)


# enh_ranked = rownames(raw)[order(rowSums(raw), decreasing = TRUE)]
# head(enh_ranked)

seurat = CreateSeuratObject(counts = raw_sparse, project = "EpiLC_enhancers")

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
nCount_violin_clusters = VlnPlot(seurat, group.by = "seurat_clusters", features = "nCount_RNA", pt.size = 0.1) +
  scale_fill_manual(values = c(brewer.pal(12, "Set3"), "darkred")) +
  ggtitle("Number of H3.3 counts over EpiLC K27ac cluster 1") +
  ylim(0, 100) +
  xlab("cluster") + 
  ylab("read count") +
  #ylim(0, 60000) +
  #scale_y_continuous(breaks = seq(0, 1000, 200)) +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 20, color = "black", angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 20, color = "black")
  ) +
  NoLegend() +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", label.y = 97)

nCount_violin_clusters

library("RColorBrewer")
brewer.pal(12, "Set3")

dim = DimPlot(object = seurat, label = FALSE, pt.size = 2, label.size = 7) + 
  scale_color_manual(values = c(brewer.pal(12, "Set3"), "darkred")) +
  xlim(-8, 8) + 
  ylim(-8, 8) + 
  ggtitle("H3.3 sciTIP-Seq UMAP - EpiLC K27ac cluster 1") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )
dim

ggarrange(dim, nCount_violin_clusters)

ggsave(
  glue("{result_folder}Seurat_H3.3_sciTIP_0510-K27ac_cl1_regions.png"),
  plot = last_plot(),
  width = 12,
  height = 5,
  dpi = 300,
)

ggsave(
  glue("{result_folder}Seurat_H3.3_sciTIP_0510-K27ac_cl1_regions.pdf"),
  plot = last_plot(),
  width = 12,
  height = 5,
  device = "pdf"
)

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
  ggtitle("read counts at EpiLC K27ac cluster 1") +
  xlab("sample") + 
  ylab("read count") + 
  #ylim(0, 100) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 14, color = "black", angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoLegend() +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", label.y = 97)
nCount_violin_samples

ggsave(
  glue("{result_folder}Seurat_H3.3_sciTIP_0510-K27ac_cl1_regions-timepoints.png"),
  plot = last_plot(),
  width = 6,
  height = 6,
  dpi = 300,
)

ggsave(
  glue("{result_folder}Seurat_H3.3_sciTIP_0510-K27ac_cl1_regions-timepoints.pdf"),
  plot = last_plot(),
  width = 6,
  height = 6,
  device = "pdf"
)

# marker analysis
marker_analysis = FindMarkers(seurat, test.use = "wilcox",
                              ident.1 = "EpiLC_6h",
                              ident.2 = "EpiLC_48h",
                              group.by = "sample",
                              min.pct = 0)

all_marker_analysis = FindAllMarkers(
  seurat,
  test.use = "roc",
  logfc.threshold = 0,
  min.pct = 0
)
