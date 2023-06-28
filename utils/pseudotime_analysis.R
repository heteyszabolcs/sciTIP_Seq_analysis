suppressPackageStartupMessages({
  library("Seurat")
  library("Signac")
  library("dyno")
  library("tidyverse")
  library("dynwrap")
  library("MatrixExtra")
})

set.seed(42)

# result folder
result_folder = "../results/Seurat/"

# load Seurat object
seurat = readRDS(file = "../data/20230316_H3.2/count_tables/20230316_H3.2_read_counts-cells_above_1000reads.Rds")

# identify top variable features
seurat = FindTopFeatures(seurat, min.cutoff = 500)
features = VariableFeatures(object = seurat)

# trajectory inference (TI)
# dynverse workflow (https://dynverse.org/users/2-quick_start/)
# counts and expression inputs: rows are cell ids, columns are bins
dataset = wrap_expression(
  counts = t_shallow(seurat@assays$sciTIP_Seq_H3.2@counts[features,]),
  expression = t_shallow(seurat@assays$sciTIP_Seq_H3.2@data[features,])
)
dataset = add_grouping(dataset,
                       seurat$seurat_clusters)
dataset = dynwrap::add_dimred(
  dataset,
  seurat@reductions$umap@cell.embeddings
)

# guidelines = guidelines_shiny(dataset)
# methods_selected = guidelines$methods_selected

# Docker must be available!
# test_docker_installation(detailed = TRUE)
# test_singularity_installation(detailed = TRUE)
model = infer_trajectory(dataset, ti_slingshot())
saveRDS(model, "/proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/results/pseudotime_analysis/slingshot_TI.Rds")

model = readRDS("../results/pseudotime_analysis/slingshot_TI.Rds")

plot_dimred(
  model, 
  expression_source = t_shallow(seurat@assays$sciTIP_Seq_H3.2@data[features,]), 
  grouping = seurat$seurat_clusters
) +
  scale_color_brewer(palette = "Set3") 

ggsave(
  "../results/pseudotime_analysis/slingshot_TI.pdf",
  plot = last_plot(),
  width = 7,
  height = 7,
  device = "pdf"
)



