library("data.table")
library("tidyverse")
library("R.utils")

# turn off warnings
options(warn = -1)

folder = "../data/coverage_win5k/"
counts = list.files(path = folder, pattern = "*.txt.gz", full.names = TRUE)

retrieve_readcount = function(count_table_path) {
  gunzip(count_table_path, remove = FALSE)
  count_table = fread(paste0(strsplit(count_table_path, split = ".txt.gz")[[1]][1], ".txt"))
  mat = as.matrix(count_table[, 5:length(count_table)])
  sums = colSums(mat)
  return(sums)
}

all_wells = unlist(lapply(counts, retrieve_readcount))
all_wells_df = tibble("overal_read_count" = all_wells, "well" = names(all_wells))

system("mkdir -p ../data/count_tables")
write_tsv(all_wells_df, "../data/count_tables/2023apr_scTIP_read_sums.tsv")

system("rm ../data/coverage_win5k/*.txt")