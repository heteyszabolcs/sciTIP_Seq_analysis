suppressPackageStartupMessages({
  library("Seurat")
  library("Signac")
  library("dyno")
  library("tidyverse")
  library("dynwrap")
  library("MatrixExtra")
  library("RColorBrewer")
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
# model = infer_trajectory(dataset, ti_slingshot())
# saveRDS(model, "/proj/snic2020-6-3/SZABOLCS/sciTIP-Seq/results/pseudotime_analysis/slingshot_TI.Rds")

model = readRDS("../results/pseudotime_analysis/slingshot_TI.Rds")

h32_pseudotime_scores = tibble(cell_id = names(calculate_pseudotime(model)), pseudotime_score = calculate_pseudotime(model))
write_tsv(h32_pseudotime_scores, "../results/pseudotime_analysis/slingshot_pseudotime_scores.tsv")

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
ggsave(
  "../results/pseudotime_analysis/slingshot_TI.png",
  plot = last_plot(),
  width = 7,
  height = 7,
  dpi = 500
)


plot_dimred(model, "pseudotime", pseudotime = calculate_pseudotime(model)) +
  ggtitle("")

ggsave(
  "../results/pseudotime_analysis/slingshot_TI_pseudotime_col.pdf",
  plot = last_plot(),
  width = 7,
  height = 7,
  device = "pdf"
)
ggsave(
  "../results/pseudotime_analysis/slingshot_TI_pseudotime_col.png",
  plot = last_plot(),
  width = 7,
  height = 7,
  dpi = 500
)

plot_dimred(model, grouping = group_onto_nearest_milestones(model)) + 
  ggtitle("Milestones")
ggsave(
  "../results/pseudotime_analysis/slingshot_TI_milestone_col.pdf",
  plot = last_plot(),
  width = 7,
  height = 7,
  device = "pdf"
)
ggsave(
  "../results/pseudotime_analysis/slingshot_TI_milestone_col.png",
  plot = last_plot(),
  width = 7,
  height = 7,
  dpi = 500
)

# early and late regions
early_h32 = read_tsv("../data/bed/early_H3.2_regions.bed", col_names = FALSE)
early_h32 = early_h32 %>% mutate(range = paste(X1, X2, X3, sep = "-")) %>% pull(range)
late_h32 = read_tsv("../data/bed/late_H3.2_regions.bed", col_names = FALSE)
late_h32 = late_h32 %>% mutate(range = paste(X1, X2, X3, sep = "-")) %>% pull(range)
mid_h32 = read_tsv("../data/bed/mid_H3.2_regions.bed", col_names = FALSE)
mid_h32 = mid_h32 %>% mutate(range = paste(X1, X2, X3, sep = "-")) %>% pull(range)

new_cols = gsub("_", "-", colnames(dataset$expression))
colnames(dataset$expression) = new_cols

FeaturePlot(object = seurat, features = sample(intersect(colnames(dataset$expression), early_h32), 9))
ggsave(
  "../results/pseudotime_analysis/FeaturePlot_earlyH3.2_ranges.pdf",
  plot = last_plot(),
  width = 12,
  height = 12,
  device = "pdf"
)
ggsave(
  "../results/pseudotime_analysis/FeaturePlot_earlyH3.2_ranges.png",
  plot = last_plot(),
  width = 12,
  height = 12,
  dpi = 500
)


FeaturePlot(object = seurat, features = sample(intersect(colnames(dataset$expression), late_h32), 9))
ggsave(
  "../results/pseudotime_analysis/FeaturePlot_lateH3.2_ranges.pdf",
  plot = last_plot(),
  width = 12,
  height = 12,
  device = "pdf"
)
ggsave(
  "../results/pseudotime_analysis/FeaturePlot_lateH3.2_ranges.png",
  plot = last_plot(),
  width = 12,
  height = 12,
  dpi = 500
)

FeaturePlot(object = seurat, features = sample(intersect(colnames(dataset$expression), mid_h32), 9))
ggsave(
  "../results/pseudotime_analysis/FeaturePlot_midH3.2_ranges.pdf",
  plot = last_plot(),
  width = 12,
  height = 12,
  device = "pdf"
)
ggsave(
  "../results/pseudotime_analysis/FeaturePlot_midH3.2_ranges.png",
  plot = last_plot(),
  width = 12,
  height = 12,
  dpi = 500
)

plot_heatmap(model, expression_source = dataset)
ggsave(
  "../results/pseudotime_analysis/slingshot_TI_heatmap.pdf",
  plot = last_plot(),
  width = 8,
  height = 6,
  device = "pdf"
)
ggsave(
  "../results/pseudotime_analysis/slingshot_TI_heatmap.png",
  plot = last_plot(),
  width = 8,
  height = 6,
  dpi = 500
)

plot_heatmap(model, expression_source = dataset$expression, 
             features_oi = sample(intersect(colnames(dataset$expression), early_h32), 15)) 
ggsave(
  "../results/pseudotime_analysis/slingshot_TI_earlyH3.2_heatmap.pdf",
  plot = last_plot(),
  width = 8,
  height = 6,
  device = "pdf"
)
ggsave(
  "../results/pseudotime_analysis/slingshot_TI_earlyH3.2_heatmap.png",
  plot = last_plot(),
  width = 8,
  height = 6,
  dpi = 500
)


plot_heatmap(model, expression_source = dataset$expression, 
             features_oi = sample(intersect(colnames(dataset$expression), late_h32), 15))
ggsave(
  "../results/pseudotime_analysis/slingshot_TI_lateH3.2_heatmap.pdf",
  plot = last_plot(),
  width = 8,
  height = 6,
  device = "pdf"
)
ggsave(
  "../results/pseudotime_analysis/slingshot_TI_lateH3.2_heatmap.png",
  plot = last_plot(),
  width = 8,
  height = 6,
  dpi = 500
)

overall_feature_importances = 
  dynfeature::calculate_overall_feature_importance(model, expression_source = dataset$expression)

features = overall_feature_importances %>% 
  top_n(40, importance) %>% 
  pull(feature_id)

branch_feature_importance = calculate_branch_feature_importance(model,
                                                                expression_source = dataset$expression)
features = branch_feature_importance %>% 
  filter(to == "M2") %>% 
  top_n(10, importance) %>% 
  pull(feature_id)

seurat = readRDS(file = "../data/20230316_H3.2/count_tables/20230316_H3.2_read_counts-cells_above_1000reads_5k.Rds")
meta = seurat@meta.data
cell_ids = rownames(meta)
meta = meta %>% mutate(cell_id = cell_ids) %>% inner_join(., h32_pseudotime_scores, by = c("cell_id" = "cell_id"))
meta = as.data.frame(meta)
rownames(meta) = cell_ids
seurat@meta.data = meta

FeaturePlot(object = seurat, features = 'pseudotime_score') +
  xlim(-10, 10) + 
  ylim(-10, 10) + 
  ggtitle("Pseudotime scores") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )

ggsave(
  glue("{result_folder}Seurat_H3.2_sciTIP_UMAP_pseudotimescores.png"),
  plot = last_plot(),
  width = 10,
  height = 10,
  dpi = 300,
)

ggsave(
  glue("{result_folder}Seurat_H3.2_sciTIP_UMAP_pseudotimescores.pdf"),
  plot = last_plot(),
  width = 10,
  height = 10,
  device = "pdf"
)
