#!/usr/bin/env Rscript
library(seqinr)
suppressMessages(library(tidyverse))

args <- commandArgs(trailingOnly = T)

fasta <- read.fasta(args[1], as.string = T, seqonly = T)

unlist(fasta) %>%
  nchar() %>%
  as_tibble() -> fasta_frame

ggplot(fasta_frame, aes(x = value)) +
  geom_histogram(
    color = "#00abff",
    binwidth = 10,
    fill = "#dbdbdb",
    alpha = 0.5
  ) +
  geom_vline(
    aes(
      xintercept = mean(value),
      color = "mean"
    ),
    linetype = "dashed",
    linewidth = 1
  ) +
  geom_vline(
    aes(
      xintercept = median(value),
      color = "median"
    ),
    linetype = "longdash",
    linewidth = 1
  ) +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_color_manual(
    name = "statistics",
    values = c(median = "#6B5CA5", mean = "#F56960")
  ) +
  ggtitle(sprintf("%s Contig Distribution", args[2])) -> plot_out

ggsave(sprintf("%s_contig_dist.png", args[2]),
  width = 8,
  height = 4,
  bg = "#FFFFFF"
)
pdf(NULL)

save(fasta_frame, file = sprintf("%s_plot_data.rdata", args[2]))
