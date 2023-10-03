library("tidyverse")
library("data.table")

# result folder
result_folder = "../results/SNPsplit/"

process_yaml = function(yaml) {
  yaml = fread(yaml, fill = TRUE)
  colnames(yaml) = "summary"
  genome1 = as.numeric(strsplit(yaml$summary[grep("^g1:", yaml$summary)], " ")[[1]][2])
  genome2 = as.numeric(strsplit(yaml$summary[grep("^g2:", yaml$summary)], " ")[[1]][2])
  unassigned = as.numeric(strsplit(yaml$summary[grep("^unassignable:", yaml$summary)], " ")[[1]][2])
  cell_id = as.character(strsplit(yaml$summary[grep("^tagged_infile:", yaml$summary)], " ")[[1]][2])
  cell_id = strsplit(cell_id, ".sam.final.allele_flagged.bam")[[1]][1]
  values = c(genome1, genome2, unassigned)
  split = c("Genome1", "Genome2", "unassigned")
  tibble = tibble(sum = values, SNPsplit = split, cell_id = cell_id)
  return(tibble)
}

total_reads = function(yaml) {
  yaml = fread(yaml, fill = TRUE)
  colnames(yaml) = "summary"
  tr = as.numeric(strsplit(yaml$summary[grep("^total_reads:", yaml$summary)], " ")[[1]][2])
  cell_id = as.character(strsplit(yaml$summary[grep("^tagged_infile:", yaml$summary)], " ")[[1]][2])
  cell_id = strsplit(cell_id, ".sam.final.allele_flagged.bam")[[1]][1]
  tibble = tibble(total_reads = tr, cell_id = cell_id)
  return(tibble)
}

yamls = list.files("../results/SNPsplit/H3.2_SNPsplit/", pattern = "yaml", full.names = TRUE)

snpsplit_output = lapply(yamls, process_yaml)
snpsplit_output = rbindlist(snpsplit_output)

abundant = snpsplit_output %>% filter(sum > 1000) %>% pull(cell_id) %>% unique %>% length
abundant

bp_input = snpsplit_output %>% 
  mutate(x_label = case_when(
    SNPsplit == "Genome1" ~ "C57BL/6J",
    SNPsplit == "Genome2" ~ "CAST/EiJ",
    .default = "Unassigned"
  ))
x_order = factor(bp_input$x_label, levels = c("Unassigned", "C57BL/6J", "CAST/EiJ"))

bp1 = ggplot(bp_input, aes(x = x_order, y = sum)) + 
  geom_boxplot(color = "red", outlier.colour = "red") +
  ylim(0, 150) +
  scale_color_grey() + 
  theme_classic() +
  labs(
    title = "SNPsplit assignments",
    x = "",
    y = "read count",
    fill = ""
  ) +
  theme(
    text = element_text(size = 13),
    plot.title = element_text(size = 17),
    axis.text.x = element_text(size = 13, color = "black", angle = 0),
    axis.title.y = element_text(size = 15, color = "black"),
    axis.text.y = element_text(size = 15, color = "black")
  ) 
bp1

ggsave(
  glue("{result_folder}SNPsplit-genome_distr_pb.pdf"),
  plot = bp1,
  width = 4,
  height = 4
)

read_counts = lapply(yamls, total_reads)
read_counts = rbindlist(read_counts)
read_counts = read_counts %>% 
  mutate(label = ifelse(total_reads == max(.$total_reads), cell_id, "")) %>% 
  mutate(feature = "total readcount")

bp = ggplot(read_counts, aes(x = feature, y = total_reads)) + 
  geom_boxplot(color = "black") +
  ylim(0, 150) +
  scale_color_grey() + 
  theme_classic() +
  labs(
    title = "",
    x = "",
    y = "",
    fill = ""
  ) +
  theme(
    text = element_text(size = 13),
    plot.title = element_text(size = 17),
    axis.text.x = element_text(size = 15, color = "black", angle = 0),
    axis.title.y = element_text(size = 15, color = "black"),
    axis.text.y = element_text(size = 15, color = "black")
  ) 
bp


yaml = "../results/SNPsplit/H3.2_SNPsplit/Cell_1A_180_S180_L001.sam.final.SNPsplit_report.yaml"
yaml = fread(yaml, fill = TRUE)
