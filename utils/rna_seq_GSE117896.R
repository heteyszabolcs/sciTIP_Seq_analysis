if (!require("pacman"))
  install.packages("pacman")
pacman::p_load("tidyverse",
               "data.table",
               "glue",
               "viridis")

# result folder
result_folder = "../results/bulk_RNA-Seq/"

# Zhifen et al. (pluripotency gene set)
# paper: https://pubmed.ncbi.nlm.nih.gov/37938437/
zhifen = fread("../data/Zhifen_et_al-pluripotency_gene_families_Fig1D.txt",
               header = FALSE)

# bulk RNA-Seq outputs (by Nextflow)
# source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117896
stringtie_outputs = list.files(
  "../data/bulk_RNA-Seq/Yang_et_al_EpiLC_GSE117896/stringtie_outputs/",
  pattern = "*abundance*",
  full.names = TRUE
)
stringtie_folder = "../data/bulk_RNA-Seq/Yang_et_al_EpiLC_GSE117896/stringtie_outputs/"

# collect expression levels
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
  Reduce(function(dtf1, dtf2)
    full_join(dtf1, dtf2, by = "gene_name"), .)
#read_counts = read_counts %>% drop_na()
read_counts = read_counts %>% group_by(gene_name) %>% summarise(across(starts_with("EpiLC"), mean)) %>%
  select(gene_name,
         EpiLC_1h,
         EpiLC_6h,
         EpiLC_12h,
         EpiLC_24h,
         EpiLC_36h,
         EpiLC_48h,
         EpiLC_72h)

write_tsv(read_counts,
          "../results/bulk_RNA-Seq/Yang_et_al-EpiLC_RNA-Seq_read_counts.tsv")

# retrieve TPMs
tpms = lapply(stringtie_outputs, collect_tpms)
tpms = tpms %>%
  Reduce(function(dtf1, dtf2)
    full_join(dtf1, dtf2, by = "gene_name"), .)
tpms = tpms %>% drop_na()
tpms = tpms %>% group_by(gene_name) %>% summarise(across(starts_with("TPM"), mean)) %>%
  select(
    gene_name,
    "TPM-EpiLC_1h",
    "TPM-EpiLC_6h",
    "TPM-EpiLC_12h",
    "TPM-EpiLC_24h",
    "TPM-EpiLC_36h",
    "TPM-EpiLC_48h",
    "TPM-EpiLC_72h"
  )

write_tsv(tpms,
          "../results/bulk_RNA-Seq/Yang_et_al-EpiLC_RNA-Seq_TPMs.tsv")

### TPM levels of a list of genes
# Yang et al. EpiLC RNA-Seq
stringtie_outputs = list.files(stringtie_folder, full.names = TRUE, pattern = "*.gene.abundance")
timepoints = c("1h", "6h", "12h", "24h", "36h", "48h", "72h")

count_tables = lapply(stringtie_outputs, function(x) {
  table = fread(x)
  return(table %>% select(`Gene ID`, TPM))
})

extract_tpm = function(gene, timepoint) {
  count_table = count_tables[grep(timepoint, x = stringtie_outputs)]
  count_table = as_tibble(count_table[[1]])
  count_table = count_table %>% filter(`Gene ID` == gene) %>% mutate(timepoint = timepoint)
  return(count_table)
}

zhifen_tpms = lapply(timepoints, function(x) {
  lapply(zhifen$V1, function(y)  {
    extract_tpm(gene = y, timepoint = x)
  })
})
zhifen_tpms = bind_rows(zhifen_tpms)
# modify gene symbols of Pou family
zhifen_tpms = zhifen_tpms %>% mutate(`Gene ID` = str_replace(zhifen_tpms$`Gene ID`, "Pou5f1", "Oct4"))
zhifen_tpms = zhifen_tpms %>% mutate(`Gene ID` = str_replace(zhifen_tpms$`Gene ID`, "Pou3f1", "Oct6"))
zhifen_tpms = zhifen_tpms %>% mutate(`Gene ID` = str_replace(zhifen_tpms$`Gene ID`, "Pou2f3", "Oct11"))

# heatmap
y_order = zhifen_tpms %>% group_by(`Gene ID`) %>% summarise(sd = sd(TPM)) %>% arrange(sd) %>% pull(`Gene ID`) %>% unique
y_order = factor(zhifen_tpms$`Gene ID`, levels = y_order)
x_order = factor(zhifen_tpms$timepoint, levels = timepoints)

hm = ggplot(zhifen_tpms, aes(x = x_order, y = y_order, fill = log2(TPM +
                                                                     1))) +
  geom_tile(color = "white",
            lwd = 0.5,
            linetype = 1) +
  scale_fill_viridis(option = "magma") +
  xlab(label = "") +
  ylab(label = "") +
  labs(fill = "log2(TPM)", title = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      color = "black",
      size = 13,
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title = element_text(size = 14)
  ) +
  coord_fixed()
print(hm)

ggsave(
  glue("{result_folder}Zhifen_et_al-pluripotency_genes_TPMS.pdf"),
  plot = hm,
  width = 5,
  height = 7
)
