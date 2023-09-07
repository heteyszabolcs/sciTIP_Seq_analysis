# packages
suppressPackageStartupMessages({
  library("Seurat")
  library("Signac")
  library("glue")
  library("tidyverse")
  library("Matrix")
  library("data.table")
  library("ggpubr")
  library("purrr")
  library("SnapATAC")
})

set.seed(42)

# result folder
result_folder = "../results/pseudotime_analysis/progression_heatmaps/"

# slingshot model
model = readRDS("../results/pseudotime_analysis/slingshot_TI.Rds")

# load Seurat object 
seurat =
  readRDS(file = "../data/20230316_H3.2/count_tables/20230316_H3.2_read_counts-cells_above_1000reads.Rds")

counts = t(seurat@assays$sciTIP_Seq_H3.2@counts)

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

# binarize H3.2 sciTIP-Seq data
snap = makeBinary(snap, mat = "bmat")

# retrieve slingshot progression
progression = model$progressions
progression = progression %>% dplyr::filter(cell_id %in% rownames(snap@bmat))

# pseudoordered cells - bins, values: snap binarized H3.2 signal
pseudoordered = snap@bmat[progression[order(progression$percentage, decreasing = FALSE),]$cell_id,
                          grepl("chr9-", colnames(snap@bmat), fixed = TRUE)]

# pseudoordered cells - bins, values: norm. H3.2 counts
pseudoordered_seurat = counts[progression[order(progression$percentage, decreasing = FALSE), ]$cell_id,
                                 grepl("chr9-", colnames(counts), fixed = TRUE)]

# create heatmap
# y axis: cells are ordered by progression along the pseudotime axis (slingshot model)
# x axis is a chromosome binned in 5 kb
library("ComplexHeatmap")
library("circlize")
col_fun = colorRamp2(c(0, 1), c("white", "red"))
parts = seq(1, dim(pseudoordered)[2], 4000)
for(i in seq(1, length(parts))) {
  start = parts[i]
  end = parts[i+1]
  pdf(
    file = glue("{result_folder}pseudotime_ordered_H3.2bin_map_{as.character(i)}_chr9.pdf"),
    width = 3,
    height = 3
  )
  print(glue("Working on: {as.character(start)}:{as.character(end)}"))
  mat = as.matrix(pseudoordered[,start:end])
  hm = Heatmap(
    mat,
    column_title = "",
    row_title = "",
    name = "replication state",
    col = col_fun,
    show_column_dend = FALSE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_row_dend = FALSE,
    #heatmap_width = unit(3, "cm"),
    #heatmap_height = unit(6, "cm"),
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_heatmap_legend = FALSE,
    use_raster = TRUE
  )
  print(hm)
  dev.off()
}