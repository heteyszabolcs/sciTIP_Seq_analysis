if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse",
               "data.table",
               "dplyr",
               "glue")

# folders
folder = "../data/20230316_H3.2/"
output = "../data/20230316_H3.2/count_tables/"

# retrieve relevant cells
read_sums = fread("../data/20230316_H3.2/count_tables/20230316_H3.2_scTIP_read_sums.tsv")
cells_above_1000reads = read_sums$well[which(read_sums$overal_read_count > 1000)]

# 80k bins
count_tables_80k = list.files(folder, pattern = ".gz", full.names = TRUE)
y = read.table(gzfile(glue("{count_tables_80k[1]}")))
ranges = y %>% mutate(range = paste(CHR, START, END, sep = "_")) %>% pull(range)

# helper function
process_counttable = function(path_to_count_table) {
  print(as.character(which(count_tables_80k == path_to_count_table)))
  x = read.table(gzfile(glue("{path_to_count_table}")))
  x = x %>% select(-CHR, -START, -END)
  return(x)
} 

# run processing and merge
merged_cts = lapply(count_tables_80k, process_counttable)
merged_cts = bind_cols(merged_cts)

# keep cells above 1000 read counts
filtered_merged_cts = merged_cts[,cells_above_1000reads]
filtered_merged_cts$range = ranges

write_tsv(filtered_merged_cts, 
          glue("{output}20230316_H3.2_read_counts-cells_above_1000reads_80k_bins.tsv"))

