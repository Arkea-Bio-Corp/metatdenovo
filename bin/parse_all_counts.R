#!/usr/bin/env Rscript
suppressMessages(library(tidyverse))
library(yaml)

args <- commandArgs(trailingOnly = TRUE)
counts_yaml <- read_yaml(args[1])

order_pipeline <- c("TRIMMOMATIC", "BOWTIE2_ALIGN", "SORTMERNA", "KRKN_NO_ARCH", 
                    "CAT_FASTQ", "BBMAP_MERGE", "BBMAP_DEDUPE")

sum_by_name <- function(count_name, yaml) {
  return(Reduce("+", yaml[grep(count_name, names(yaml))]))
}

res <- data.frame(tool   = order_pipeline,
                  counts = c(sum_by_name("TRIMMOMATIC", counts_yaml),
                             sum_by_name("BOWTIE", counts_yaml),
                             sum_by_name("SORTMERNA", counts_yaml),
                             sum_by_name("KRKN_NO_ARCH", counts_yaml),
                             counts_yaml$`NFCORE_METATDENOVO:METATDENOVO:CAT_FASTQ`,
                             counts_yaml$`NFCORE_METATDENOVO:METATDENOVO:BBMAP_MERGE`,
                             counts_yaml$`NFCORE_METATDENOVO:METATDENOVO:BBMAP_DEDUPE`),
                  effect = c("Trim", "Filter", "Filter", "Filter", "Combine", "Merge", 
                             "Filter"))
res$tool <- factor(res$tool, levels=order_pipeline)

ggplot(res, aes(x=tool, y=counts, fill = effect)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  ggtitle("Total Counts After Each Step, before Assembly") +
  xlab("Tool") +
  ylab("Number of reads") +
  scale_y_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6)) +
  scale_fill_manual("Effect", values = c("#E03616", "#5C946E", "#F45B69", "#645E9D")) +
  theme(axis.text.x = element_text(angle = 35, vjust = 0.75),
        legend.box.background = element_rect(colour = "black"),
        legend.position = "top")

ggsave(filename = "All_counts_plot_mqc.png", device = "png",
       width = 8, height = 4, units = "in", bg = "#FFFFFF")
