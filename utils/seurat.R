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
  library("SnapATAC")
})

set.seed(42)

# EpiLC enhancer related cell ids
ids = fread("../results/Seurat/cells_with_H3.3_over_EpiLC_enhancers.tsv")
ids2 = fread("../results/Seurat/cells_with_highH3.3_over_EpiLC_enhancers.tsv")

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
# 
# # export Rds
# saveRDS(seurat, "../data/20230316_H3.2/count_tables/20230316_H3.2_read_counts-cells_above_1000reads.Rds")

# load Seurat object
seurat = readRDS(file = "../data/20230510_EpiLC/count_tables/20230510_EpiLC_read_counts-cells_above_1000reads.Rds")

# normalization
seurat = RunTFIDF(seurat)

# dubstepR.out = DUBStepR(input.data = seurat@assays$RNA@data, 
#                         min.cells = 0.05*ncol(seurat), 
#                         optimise.features = T, 
#                         k = 10, 
#                         species = "mouse",
#                         num.pcs = 20, error = 0)
# seurat@assays$RNA@var.features = dubstepR.out$optimal.feature.genes

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

FeaturePlot(seurat, 
            features = dubstepR.out$optimal.feature.genes, 
            cols = c("lightgrey", "magenta"))

# quality plots
nCount_violin_clusters = VlnPlot(seurat, group.by = "seurat_clusters", features = "nCount_RNA", pt.size = 0.1) +
  scale_fill_brewer(palette = "Set3") +
  ggtitle("nCount (bins)") +
  xlab("cluster") + 
  ylab("read count") +
  #ylim(0, 60000) +
  scale_y_continuous(breaks = seq(0, 60000, 10000)) +
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
  ggtitle("H3.3 sciTIP-Seq") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )
dim

ggsave(
  glue("{result_folder}Seurat_H3.3_sciTIP_0510_UMAP.png"),
  plot = dim,
  width = 10,
  height = 10,
  dpi = 300,
)

ggsave(
  glue("{result_folder}Seurat_H3.3_sciTIP_0510_UMAP.pdf"),
  plot = dim,
  width = 10,
  height = 10,
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
  ggtitle("nCount (bins)") +
  xlab("sample") + 
  ylab("read count") + 
  #ylim(0, 60000) +
  scale_y_continuous(breaks = seq(0, 60000, 10000)) +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 14, color = "black", angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 25, color = "black")
  ) +
  NoLegend()
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

ggsave(
  glue("{result_folder}Seurat_H3.3_sciTIP_0510-sample_level_UMAP.png"),
  plot = dim_samples,
  width = 10,
  height = 10,
  dpi = 300,
)

ggsave(
  glue("{result_folder}Seurat_H3.3_sciTIP_0510-sample_level_UMAP.pdf"),
  plot = dim_samples,
  width = 10,
  height = 10,
  device = "pdf"
)

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
  glue("{result_folder}Seurat_H3.3_sciTIP_0510-quality_plots.png"),
  plot = violins,
  width = 12,
  height = 5,
  dpi = 300,
)

ggsave(
  glue("{result_folder}Seurat_H3.3_sciTIP_0510-quality_plots.pdf"),
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
                              logfc.threshold = 0.1)
summary = marker_analysis %>% 
  dplyr::filter(abs(avg_log2FC) > 0.1) %>% 
  dplyr::filter(p_val_adj < 0.05)

# on whole dataset
all_marker_analysis = FindAllMarkers(
  seurat,
  test.use = "roc",
  group.by = "sample",
  logfc.threshold = 0.1
)

# cells with H3.3 occupied EpiLC spec. enhancers
ids = ids$cells_with_H3.3_over_EpiLC_enhancers
ids = paste0(ids, "_L001_R1_001.fastq")

DimPlot(
  object = seurat,
  pt.size = 2,
  cells.highlight = ids,
  cols.highlight = "red",
  cols = "gray",
  order = TRUE
) +
  xlim(-10, 10) +
  ylim(-10, 10) +
  scale_color_manual(
    values = c("#bdbdbd", "#de2d26"),
    labels = c("cells with H3.3-less EpiLC enhancers", "cells with H3.3 EpiLC enhancers")
  ) 

ggsave(
  glue("{result_folder}cells_w_H3.3EpiLC_enh-UMAP.png"),
  plot = last_plot(),
  width = 10,
  height = 10,
  dpi = 300
)

ggsave(
  glue("{result_folder}cells_w_H3.3EpiLC_enh-UMAP.pdf"),
  plot = last_plot(),
  width = 10,
  height = 10,
  device = "pdf"
)

# cells with high H3.3 occupied EpiLC spec. enhancers
ids2 = ids2$cells_with_highH3.3_over_EpiLC_enhancers
ids2 = paste0(ids2, "_L001_R1_001.fastq")

DimPlot(
  object = seurat,
  pt.size = 2,
  cells.highlight = ids2,
  cols.highlight = "red",
  cols = "gray",
  order = TRUE
) +
  xlim(-10, 10) +
  ylim(-10, 10) +
  scale_color_manual(
    values = c("#bdbdbd", "#de2d26"),
    labels = c("cells with low/zero H3.3 EpiLC enhancers", "cells with high H3.3 EpiLC enhancers")
  ) 

ggsave(
  glue("{result_folder}cells_w_highH3.3EpiLC_enh-UMAP.png"),
  plot = last_plot(),
  width = 10,
  height = 10,
  dpi = 300
)

ggsave(
  glue("{result_folder}cells_w_highH3.3EpiLC_enh-UMAP.pdf"),
  plot = last_plot(),
  width = 10,
  height = 10,
  device = "pdf"
)

# module scores based on k27ac clusters of EpiLC clusters
# https://satijalab.org/seurat/reference/addmodulescore
# source: Tirosh et al, Science (2016)
# the module score represents relative expression (relative H3.3 occupancy). If bin A is highly expressed across all cells, 
# the module scores assess whether a given cell expresses bin A more often than other highly expressed bins 
# Therefore, if the module score is high, it indicates these bins have high h3.3. 
# If the module score is low, it might mean that these bins have h3.3 but not more so than in other cells.
k27ac_cluster2 = fread("../data/bed/Yang_EpiLC_K27ac_over_epilc_enh-kmeans-cluster2.bed")
k27ac_cluster2 = GRanges(
  seqnames = k27ac_cluster2$V1,
  ranges = IRanges(
    start = k27ac_cluster2$V2,
    end = k27ac_cluster2$V3
  )
)
k27ac_cluster1 = fread("../data/bed/Yang_EpiLC_K27ac_over_epilc_enh-kmeans-cluster1.bed")
k27ac_cluster1 = GRanges(
  seqnames = k27ac_cluster1$V1,
  ranges = IRanges(
    start = k27ac_cluster1$V2,
    end = k27ac_cluster1$V3
  )
)

bins = tibble(regions = rownames(seurat@assays$RNA@data)) %>% separate(., regions, sep = "-", into = c("V1", "V2", "V3")) %>% 
  mutate(V2 = as.numeric(V2), V3 = as.numeric(V3))
bins = GRanges(
  seqnames = bins$V1,
  ranges = IRanges(
    start = bins$V2,
    end = bins$V3
  )
)

ol2 = findOverlaps(k27ac_cluster2, bins, type = "any", ignore.strand = FALSE)
# intersection
ol2 = bins[subjectHits(ol2)]
ol2 = as_tibble(ol2)
ol2 = ol2 %>% mutate(id = paste(seqnames, start, end, sep = "-")) %>% pull(id)

seurat = AddModuleScore(
  object = seurat,
  features = list(sample(ol2, 500, replace = FALSE)),
  name = "cluster2_modulescore",
  seed = 1
)
meta = seurat@meta.data
# seurat@meta.data = meta_modscore

FeaturePlot(object = seurat, features = "cluster2_modulescore1", pt.size = 2, order = TRUE) +
  xlim(-10, 10) + 
  ylim(-10, 10) + 
  labs(title = "late H3K27ac-ed enhancers", subtitle = "k-means cluster 2") +
  guides(fill = guide_legend(title="module score")) +
  scale_colour_gradient2(low = "#f0f0f0",
                         mid = "#ffeda0",
                         high = "#de2d26",
                         midpoint = 0.075) +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20, hjust = 0),
    plot.title.position = "plot",
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )

ggsave(
  glue("{result_folder}Seurat_H3.3-enhancer_module_scores-k27ac_cl2.png"),
  plot = last_plot(),
  width = 8,
  height = 8,
  dpi = 300,
)

ggsave(
  glue("{result_folder}Seurat_H3.3-enhancer_module_scores-k27ac_cl2.pdf"),
  plot = last_plot(),
  width = 8,
  height = 8,
  device = "pdf"
)


ol1 = findOverlaps(k27ac_cluster1, bins, type = "any", ignore.strand = FALSE)
# intersection
ol1 = bins[subjectHits(ol1)]
ol1 = as_tibble(ol1)
ol1 = ol1 %>% mutate(id = paste(seqnames, start, end, sep = "-")) %>% pull(id)

seurat = AddModuleScore(
  object = seurat,
  features = ol1,
  name = "cluster1_modulescore",
  seed = 1
)
meta = seurat@meta.data
# seurat@meta.data = meta_modscore
FeaturePlot(object = seurat, features = "cluster1_modulescore1", pt.size = 2, order = TRUE) +
  xlim(-10, 10) + 
  ylim(-10, 10) + 
  labs(title = "early H3K27ac-ed enhancers", subtitle = "k-means cluster 1") +
  guides(fill = guide_legend(title="module score")) +
  scale_colour_gradient2(low = "#f0f0f0",
                         mid = "#ffeda0",
                         high = "#de2d26",
                         midpoint = 0.075) +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20, hjust = 0),
    plot.title.position = "plot",
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )

ggsave(
  glue("{result_folder}Seurat_H3.3-enhancer_module_scores-k27ac_cl1.png"),
  plot = last_plot(),
  width = 8,
  height = 8,
  dpi = 300,
)

ggsave(
  glue("{result_folder}Seurat_H3.3-enhancer_module_scores-k27ac_cl1.pdf"),
  plot = last_plot(),
  width = 8,
  height = 8,
  device = "pdf"
)

object = seurat
features = list(ol2)
pool = rownames(object)
nbin = 24
ctrl = 100
k = FALSE
name = "cluster2"
seed = 1
cluster.length <- length(x = features)
assay.data <- GetAssayData(object = object)
# For all genes, get the average expression across all cells (named vector)
data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
# Order genes from lowest average expression to highest average expression
data.avg <- data.avg[order(data.avg)]
# Use ggplot2's cut_number function to make n groups with (approximately) equal numbers of observations. 
# The 'rnorm(n = length(data.avg))/1e+30' part adds a tiny bit of noise to the data, presumably to break ties.
data.cut <- ggplot2::cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30,
                                n = nbin,
                                labels = FALSE,
                                right = FALSE)
# Set the names of the cuts as the gene names
names(x = data.cut) <- names(x = data.avg)
ctrl.use <- vector(mode = "list", length = cluster.length)
for (i in 1:cluster.length) {
  # Get the gene names from the input gene set as a character vector  
  features.use <- features[[i]]
  
  # Loop through the provided genes (1:num_genes) and for each gene, find ctrl (default=100) genes from the same expression bin (by looking in data.cut):
  for (j in 1:length(x = features.use)) {
    # Within this loop, 'data.cut[features.use[j]]' gives us the expression bin number. We then sample `ctrl` genes from that bin without replacement and add the gene names to ctrl.use.
    ctrl.use[[i]] <- c(ctrl.use[[i]],
                       names(x = sample(x = data.cut[which(x = data.cut == data.cut[features.use[j]])],
                                        size = ctrl,
                                        replace = FALSE)))
  }
}

pdf(
  file = glue("{result_folder}H3.3_enhancer_module_score_plot-k27ac_cluster2.pdf"),
  width = 6,
  height = 6
)
plot(data.avg, pch=16, ylab="Average H3.3 level across all cells", xlab = "All bins, ranked")
# Add red points for selected control genes
points(which(names(data.avg)%in%ctrl.use[[1]]), data.avg[which(names(data.avg)%in%ctrl.use[[1]])], pch=16, col="#7fcdbb")

# Add blue points for genes in the input gene list
points(which(names(data.avg)%in%features[[1]]), data.avg[which(names(data.avg)%in%features[[1]])], pch=16, col="#de2d26")

# Add a legend
legend(x = "topleft",
       legend = c("bin", "selected control bin", "bins overlap with K27ac cluster 2"),
       col = c("black", "#de2d26", "#7fcdbb"),
       pch = 16)
dev.off()

features = list(ol1)
name = "cluster1"
seed = 1
cluster.length <- length(x = features)
assay.data <- GetAssayData(object = object)
# For all genes, get the average expression across all cells (named vector)
data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
# Order genes from lowest average expression to highest average expression
data.avg <- data.avg[order(data.avg)]
# Use ggplot2's cut_number function to make n groups with (approximately) equal numbers of observations. 
# The 'rnorm(n = length(data.avg))/1e+30' part adds a tiny bit of noise to the data, presumably to break ties.
data.cut <- ggplot2::cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30,
                                n = nbin,
                                labels = FALSE,
                                right = FALSE)
# Set the names of the cuts as the gene names
names(x = data.cut) <- names(x = data.avg)
ctrl.use <- vector(mode = "list", length = cluster.length)
for (i in 1:cluster.length) {
  # Get the gene names from the input gene set as a character vector  
  features.use <- features[[i]]
  
  # Loop through the provided genes (1:num_genes) and for each gene, find ctrl (default=100) genes from the same expression bin (by looking in data.cut):
  for (j in 1:length(x = features.use)) {
    # Within this loop, 'data.cut[features.use[j]]' gives us the expression bin number. We then sample `ctrl` genes from that bin without replacement and add the gene names to ctrl.use.
    ctrl.use[[i]] <- c(ctrl.use[[i]],
                       names(x = sample(x = data.cut[which(x = data.cut == data.cut[features.use[j]])],
                                        size = ctrl,
                                        replace = FALSE)))
  }
}

pdf(
  file = glue("{result_folder}H3.3_enhancer_module_score_plot-k27ac_cluster1.pdf"),
  width = 6,
  height = 6
)
plot(data.avg, pch=16, ylab="Average H3.3 level across all cells", xlab = "All bins, ranked")
# Add red points for selected control genes
points(which(names(data.avg)%in%ctrl.use[[1]]), data.avg[which(names(data.avg)%in%ctrl.use[[1]])], pch=16, col="#7fcdbb")

# Add blue points for genes in the input gene list
points(which(names(data.avg)%in%features[[1]]), data.avg[which(names(data.avg)%in%features[[1]])], pch=16, col="#de2d26")

# Add a legend
legend(x = "topleft",
       legend = c("bin", "selected control bin", "bins overlap with K27ac cluster 1"),
       col = c("black", "#de2d26", "#7fcdbb"),
       pch = 16)
dev.off()