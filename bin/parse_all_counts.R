#!/usr/bin/env Rscript
suppressMessages(library(tidyverse))
args <- commandArgs(trailingOnly = TRUE)

counts_input <- as_tibble(read.csv(args[1],
  header = F,
  col.names = c("Tool", "Reads")
))

order_pipeline <- c(
  "TRIMMOMATIC", "BOWTIE2_ALIGN", "SORTMERNA", "KRKN_NO_ARCH",
  "CAT_FASTQ", "BBMAP_MERGE", "BBMAP_DEDUPE"
)

counts_input %>%
  separate(Tool, into = c(NA, NA, "Tool"), sep = ":", remove = FALSE) %>%
  group_by(Tool) %>%
  summarise(across(everything(), sum)) %>%
  within(Tool <- factor(Tool, levels = order_pipeline)) -> counts_table

ggplot(counts_table, aes(x = Tool, y = Reads, fill = Reads)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "darkgreen", high = "red") +
  theme_minimal() +
  ggtitle("Total Counts After Each Step, before Assembly") +
  xlab("Tool") +
  ylab("Number of reads") +
  scale_y_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6)) +
  theme(
    axis.text.x = element_text(angle = 35, vjust = 0.75),
    legend.position = "none"
  )

ggsave(
  filename = "All_counts_plot_mqc.png", device = "png",
  width = 8, height = 4, units = "in", bg = "#FFFFFF"
)
