library("data.table")
library("tidyverse")
library("R.utils") # for unzipping

# turn off warnings
options(warn = -1)

# read count table chunks (5 kb bins)
folder = "../data/20230316_H3.2/coverage_win5k/"
counts = list.files(path = folder, pattern = "*.txt.gz", full.names = TRUE)

# retrieve total read count per cell
retrieve_readcount = function(count_table_path) {
  gunzip(count_table_path, remove = FALSE)
  count_table = fread(paste0(strsplit(count_table_path, split = ".txt.gz")[[1]][1], ".txt"))
  mat = as.matrix(count_table[, 5:length(count_table)])
  sums = colSums(mat)
  return(sums)
}

all_wells = unlist(lapply(counts, retrieve_readcount))
all_wells_df = tibble("overal_read_count" = all_wells, "well" = names(all_wells))

system("mkdir -p ../data/20230316_H3.2/count_tables")
write_tsv(all_wells_df, "../data/20230316_H3.2/count_tables/20230316_H3.2_scTIP_read_sums.tsv")

system("rm ../data/20230316_H3.2/coverage_win5k/*.txt")

# retrieve mean read count per cell
retrieve_mean = function(count_table_path) {
  gunzip(count_table_path, remove = FALSE)
  count_table = fread(paste0(strsplit(count_table_path, split = ".txt.gz")[[1]][1], ".txt"))
  mat = as.matrix(count_table[, 5:length(count_table)])
  means = colMeans(mat)
  return(means)
}

all_wells = unlist(lapply(counts, retrieve_mean))
all_wells_df = tibble("mean_read_count" = all_wells, "well" = names(all_wells))

# compute whole read count table with cells having more than 1000 reads
read_tables = function(count_table_path) {
  #gunzip(count_table_path, remove = FALSE)
  count_table = fread(paste0(strsplit(count_table_path, split = ".txt.gz")[[1]][1], ".txt"))
  rows = count_table %>% mutate(rowname = paste(CHR, START, END, sep = "_")) %>% pull(rowname)
  mat = as.matrix(count_table[, 5:length(count_table)])
  mat = mat[, colSums(mat) > 1000]
  df = as.data.frame(mat)
  rownames(df) = rows
  return(df)
}

all_counts = lapply(counts, read_tables)
rows = rownames(all_counts[[1]][1])
all_counts = bind_cols(all_counts)
all_counts = all_counts %>% mutate(range = rows)


# write out for Seurat analysis
write_tsv(all_counts, "../data/20230316_H3.2/count_tables/20230316_H3.2_read_counts-cells_above_1000reads.tsv")



