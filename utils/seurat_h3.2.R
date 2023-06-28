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
raw = fread(file = "../data/20230316_H3.2/count_tables/20230316_H3.2_read_counts-cells_above_1000reads.tsv")
rows = raw$range
raw = raw %>% dplyr::select(-range)
raw_sparse = as(as.matrix(raw), "sparseMatrix")
rownames(raw_sparse) = rows
rm(raw)
seurat = CreateSeuratObject(counts = raw_sparse, project = "sciTIP_Seq")
seurat = RenameAssays(object = seurat, RNA = 'sciTIP_Seq_H3.2')

# export Rds
saveRDS(seurat, "../data/20230316_H3.2/count_tables/20230316_H3.2_read_counts-cells_above_1000reads.Rds")

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

seurat = RenameAssays(object = seurat, RNA = 'sciTIP_Seq_H3.2')
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

# coembed with H3.3
h33 = readRDS(file = "../data/count_tables/20230510_EpiLC_read_counts-cells_above_1000reads.Rds")
h33@meta.data$data_type = "H3.3"
seurat@meta.data$data_type = "H3.2"

# find top variable features for anchor finding
vf = FindTopFeatures(seurat, min.cutoff = 1000)
vf = VariableFeatures(vf)

# find anchors between to two sets
transfer.anchors = FindTransferAnchors(
  reference = h33, 
  features = vf,
  query = seurat, reduction = "cca"
)

refdata = GetAssayData(h33, assay = "sciTIP_Seq_H3.2", slot = "data")[vf, ]
imputation = TransferData(anchorset = transfer.anchors, refdata = refdata, 
                           weight.reduction = seurat[["lsi"]],  dims = 2:30)

seurat[["sciTIP_Seq_H3.2"]] = imputation
coembed = merge(x = seurat, y = h33)

# normalization
coembed = RunTFIDF(coembed)
coembed = FindTopFeatures(coembed, min.cutoff = 'q0')
coembed = RunSVD(coembed)

# Non-linear dimension reduction and clustering
coembed = RunUMAP(object = coembed,
                  reduction = 'lsi',
                  dims = 2:30)
coembed = FindNeighbors(object = coembed,
                        reduction = 'lsi',
                        dims = 2:30)
coembed = FindClusters(object = coembed,
                       verbose = FALSE,
                       resolution = 0.5,
                       algorithm = 3)

dim = DimPlot(object = coembed, group.by = "data_type", label = FALSE, pt.size = 2, label.size = 7) + 
  scale_color_manual(values = c("#fc9272", "#9ecae1")) +
  xlim(-20, 20) + 
  ylim(-20, 20) + 
  ggtitle("H3.3/H3.2 coembedding") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )
dim

# quality plots
nCount_violin_clusters = VlnPlot(coembed, group.by = "data_type", features = "nCount_sciTIP_Seq_H3.2", pt.size = 0.1) +
  scale_color_manual(values = c("#fc9272", "#9ecae1")) +
  ggtitle("") +
  xlab("") + 
  ylab("read count") +
  ylim(0, 500000) +
  #scale_y_continuous(breaks = seq(0, 500000, 100000)) +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black", angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoLegend() +
  stat_compare_means(label.y = 400000, label.x = 1.25) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "H3.3", label.y = 350000)

nCount_violin_clusters

ggsave(
  glue("{result_folder}H3.3_H3.2_coembed-readcounts.pdf"),
  plot = nCount_violin_clusters,
  width = 7,
  height = 5,
  device = "pdf"
)

ggsave(
  glue("{result_folder}H3.3_H3.2_coembed-readcounts.png"),
  plot = nCount_violin_clusters,
  width = 7,
  height = 5,
  dpi = 500
)


ggsave(
  glue("{result_folder}H3.3_H3.2_coembed-UMAP.pdf"),
  plot = dim,
  width = 7,
  height = 5,
  device = "pdf"
)

ggsave(
  glue("{result_folder}H3.3_H3.2_coembed-UMAP.png"),
  plot = dim,
  width = 7,
  height = 5,
  dpi = 500
)

