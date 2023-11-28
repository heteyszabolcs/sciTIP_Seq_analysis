if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse",
               "data.table",
               "glue") 

stringtie_outputs = list.files("../data/bulk_RNA-Seq/Yang_et_al_EpiLC_GSE117896/stringtie_outputs/",
                               pattern = "*abundance*", full.names = TRUE)
stringtie_folder = "../data/bulk_RNA-Seq/Yang_et_al_EpiLC_GSE117896/stringtie_outputs/"

collect_readcounts = function(stringtie_output) {
  colname = strsplit(stringtie_output, stringtie_folder)[[1]][2]
  colname = strsplit(colname, ".gene.abundance.txt")[[1]][1]
  print(colname)
  output = fread(stringtie_output)
  output = output %>% select(gene_name = `Gene Name`, Coverage) %>%
    mutate(Coverage = round(Coverage, 2))
  colnames(output) = c("gene_name", colname)
  return(output)
}

collect_tpms = function(stringtie_output) {
  colname = strsplit(stringtie_output, stringtie_folder)[[1]][2]
  colname = strsplit(colname, ".gene.abundance.txt")[[1]][1]
  print(colname)
  output = fread(stringtie_output)
  output = output %>% select(gene_name = `Gene Name`, TPM) %>%
    mutate(TPM = round(TPM, 2))
  colnames(output) = c("gene_name", glue("TPM-{colname}"))
  return(output)
}

# retrieve read counts
read_counts = lapply(stringtie_outputs, collect_readcounts)
read_counts = read_counts %>%
  Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="gene_name"), .)
#read_counts = read_counts %>% drop_na()
read_counts = read_counts %>% group_by(gene_name) %>% summarise(across(starts_with("EpiLC"), mean)) %>% 
  select(gene_name, EpiLC_1h, EpiLC_6h, EpiLC_12h, EpiLC_24h, EpiLC_36h, EpiLC_48h, EpiLC_72h)

write_tsv(read_counts, "../results/bulk_RNA-Seq/Yang_et_al-EpiLC_RNA-Seq_read_counts.tsv")

# retrieve read counts
tpms = lapply(stringtie_outputs, collect_tpms)
tpms = tpms %>%
  Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="gene_name"), .)
tpms = tpms %>% drop_na()
tpms = tpms %>% group_by(gene_name) %>% summarise(across(starts_with("TPM"), mean)) %>% 
  select(gene_name, "TPM-EpiLC_1h", "TPM-EpiLC_6h", "TPM-EpiLC_12h", "TPM-EpiLC_24h", "TPM-EpiLC_36h", 
         "TPM-EpiLC_48h", "TPM-EpiLC_72h")

write_tsv(tpms, "../results/bulk_RNA-Seq/Yang_et_al-EpiLC_RNA-Seq_TPMs.tsv")


