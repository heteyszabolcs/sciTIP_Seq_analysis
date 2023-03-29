library("data.table"); library("tidyverse")

allcombs = fread("../data/all_combinations.tsv")
r5 = fread("../data/r5_barcodes.txt")
i5 = fread("../data/i5_barcodes.txt")
i7 = fread("../data/i7_barcodes.txt")

find_combination = function(combindex, header = x) {
  if (str_detect(header, combindex)) {
    return(combindex)
  }
}

unlist(lapply(r5$r5_sequence, find_combination, header = "TACTTGAA+TAGGCAGCGCGTGGATCTCCCTATGCAGT"))
unlist(lapply(i5$i5_sequence, find_combination, header = "TACTTGAA+TAGGCAGCGCGTGGATCTCCCTATGCAGT"))
unlist(lapply(i7$i7_sequence, find_combination, header = "TACTTGAA+TAGGCAGCGCGTGGATCTCCCTATGCAGT"))


