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
})

set.seed(42)

# result folder
result_folder = "../results/Seurat/"

# create Seurat object
# raw = fread(file = "../data/20230316_H3.2/count_tables/20230316_H3.2_read_counts-cells_above_1000reads.tsv")
# rows = raw$range
# raw = raw %>% dplyr::select(-range)
# raw_sparse = as(as.matrix(raw), "sparseMatrix")
# rownames(raw_sparse) = rows
# rm(raw)
# seurat = CreateSeuratObject(counts = raw_sparse, project = "sciTIP_Seq")
# seurat = RenameAssays(object = seurat, RNA = 'sciTIP_Seq_H3.2')
# 
# # export Rds
# saveRDS(seurat, "../data/20230316_H3.2/count_tables/20230316_H3.2_read_counts-cells_above_1000reads.Rds")

# load Seurat object
seurat = readRDS(file = "../data/20230316_H3.2/count_tables/20230316_H3.2_read_counts-cells_above_1000reads.Rds")

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

# resolution 0.5 with original Louvain algorithm
seurat = FindClusters(object = seurat,
                         verbose = FALSE,
                         resolution = 0.5,
                         algorithm = 1)

saveRDS(seurat, file = "../data/20230316_H3.2/count_tables/20230316_H3.2_read_counts-cells_above_1000reads.Rds")

# quality plots
nCount_violin_clusters = VlnPlot(seurat, group.by = "seurat_clusters", features = "nCount_RNA", pt.size = 0.1) +
  scale_fill_brewer(palette = "Set3") +
  ggtitle("nCount (bins)") +
  xlab("cluster") + 
  ylab("read count") +
  ylim(0, 300000) +
  scale_y_continuous(breaks = seq(0, 300000, 100000)) +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black", angle = 0),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoLegend()
nCount_violin_clusters

ggsave(
  glue("{result_folder}Seurat_H3.2_sciTIP_0316-quality_plots.png"),
  plot = last_plot(),
  width = 10,
  height = 10,
  dpi = 300,
)

ggsave(
  glue("{result_folder}Seurat_H3.2_sciTIP_0316-quality_plots.pdf"),
  plot = last_plot(),
  width = 10,
  height = 10,
  device = "pdf"
)
# visualizations
dim = DimPlot(object = seurat, label = FALSE, pt.size = 2, label.size = 7) + 
  scale_color_brewer(palette = "Set3") +
  xlim(-10, 10) + 
  ylim(-10, 10) + 
  ggtitle("H3.2 sciTIP-Seq") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )
dim

meta = seurat@meta.data
for(cluster in unique(meta$seurat_clusters)) {
  df = meta %>% dplyr::filter(seurat_clusters == cluster)
  df = df %>% mutate(id = rownames(df)) %>% separate(id, sep = ".fastq", into = "id") %>% 
    dplyr::select(id)
  df = as_tibble(df)
  write_tsv(df, glue("../results/Seurat/H3.2_res0.5_cluster{as.character(cluster)}_ids.tsv"),
            col_names = FALSE)
}

ggsave(
  glue("{result_folder}Seurat_H3.2_sciTIP_0316_UMAP.png"),
  plot = dim,
  width = 10,
  height = 10,
  dpi = 300,
)

ggsave(
  glue("{result_folder}Seurat_H3.2_sciTIP_0316_UMAP.pdf"),
  plot = dim,
  width = 10,
  height = 10,
  device = "pdf"
)

# find H3.2 marker regions of the clusters
marker_analysis = FindAllMarkers(seurat)
summary = marker_analysis %>% 
  dplyr::filter(avg_log2FC > 0.2) %>% 
  dplyr::filter(p_val_adj < 0.05) %>% 
  group_by(cluster) %>% count()

write_tsv(marker_analysis, glue("{result_folder}Seurat_H3.2_sciTIP_Marker_analysis.tsv"))

# example
# cluster 1 - H3.2 marker region
FeaturePlot(seurat, features = "chr2-181930000-181935000")



