library("tidyverse")
library("data.table")
library("ggplot2")
library("ggpubr")
library("glue")

result_folder = "../results/explorative/"

cov = fread("../data/coverage_win5k/1A_sciTIP_counts_matrix.txt")
cov[1:10, 1:10]

return_max = function(x) {
  max(cov[,..x])
}
max_values = tibble(max = sapply(5:length(colnames(cov)), return_max))

max_hist = ggplot(max_values,
              aes(x = log2(max))) +
  geom_histogram(position = "identity", alpha = 0.8, fill = "#fc9272", color = "black") +
  geom_density(alpha = 1.0) +
  labs(title = "",
       x = "log2(max read count) of A1",
       y = "Density") +
  theme_classic() +
  theme(
    text = element_text(size = 18),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 15, color = "black")
  )
max_hist

mat = cov[,5:length(colnames(cov))]
mat = as.matrix(mat)
rc_per_cell = tibble(read_count = colSums(mat), well = "A1") 

violin = ggplot(rc_per_cell, aes(x = well, y = log10(read_count))) +
  geom_violin(color = "black", fill = "#fc9272") +
  scale_fill_brewer(palette = "Reds") +
  #ylim(0, 10000) +
  labs(
    title = "",
    x = "",
    y = "log10(read count)"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 10),
    axis.text.x = element_text(size = 15, color = "black"),
    axis.text.y = element_text(size = 15, color = "black")
  )
violin

ggarrange(max_hist, violin)

ggsave(
  glue("{result_folder}Well_1A_readcounts.pdf"),
  width = 12,
  height = 6
)

ggsave(
  glue("{result_folder}Well_1A_readcounts.png"),
  width = 12,
  height = 6,
  dpi = 500
)



