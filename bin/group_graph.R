library(tidyverse)
library(patchwork)
# setwd("~/Documents/metatdenovo/localdata/assembly_compare/small/")
# setwd("~/Documents/metatdenovo/localdata/assembly_compare/medium/")
# setwd("~/Documents/metatdenovo/localdata/assembly_compare/med_plass_ta/")
setwd("~/Documents/metatdenovo/localdata/assembly_compare/large/")

temp = list.files(pattern="assemblies.csv", recursive = T)
myfiles = lapply(temp, read.delim)

temp %>%
  map_df(~read_csv(.)) %>%
  cbind(temp, .) -> transrate_tab

transrate_tab %>%
  separate(1, c("Assembler"), sep = "/") -> transrate_tab

list.files(pattern = ".bowtie2.log", recursive = T) %>%
  lapply(read_lines, skip = 5) %>%
  unlist() %>%
  data.frame() %>%
  separate(1, c("."), sep = "%") -> bt2_percent

cbind(transrate_tab, sapply(bt2_percent, as.numeric)) -> stats_table
rename(stats_table, bt2 = .) %>%
  select(Assembler, bt2, n_seqs, mean_len, n50) -> stats_table

# Plot Transrate Stats ----------------------------------------------------


ggplot(stats_table, aes(Assembler, bt2)) +
  geom_bar(aes(fill = Assembler), stat = "identity") +
  scale_fill_manual(values=c("#69b3a2", "#404080", "#938BA1",
                             "#F5D491", "#B6244F", "#FF8552")) +
  lims(y = c(0, 100)) +
  geom_text(aes(y = bt2, label = sprintf("%s%%", bt2)),
            position = position_dodge(width = 0.9),
            vjust=-0.25) +
  theme_bw() +
  ggtitle("Bowtie Read Alignment Percentage") +
  theme(legend.position = "none") -> bowtie

ggplot(stats_table, aes(Assembler)) +
  geom_bar(aes(y = n_seqs, fill = Assembler), stat = "identity") +
  geom_point(shape = 17, size = 4,
             aes(y = mean_len, color = Assembler),
             stat = "identity") +
  scale_fill_manual(values=c("#69b3a2", "#404080", "#938BA1",
                             "#F5D491", "#B6244F", "#FF8552")) +
  scale_y_continuous(name = "Number of Contigs",
                     # sec.axis = dup_axis(name = "Mean Length (▲)")
                     sec.axis = sec_axis(log10, name = "Mean Length (▲)")
                     ) +
  theme_bw() +
  theme(legend.position = "bottom") +
  ggtitle("Assembler contig count and mean length") +
  geom_text(aes(y= n_seqs, label=n_seqs),
            position=position_dodge(width=0.9),
            vjust=-0.25) -> contigs

ggplot(stats_table, aes(Assembler)) +
  geom_bar(aes(y = n50, fill = Assembler), stat = "identity") +
  scale_fill_manual(values=c("#69b3a2", "#404080", "#938BA1",
                             "#F5D491", "#B6244F", "#FF8552")) +
  ggtitle("N50 value by Assembler") +
  theme_bw() +
  theme(legend.position = "none") -> n50

# Compare contig distribution ---------------------------------------------

load("Megahit/Megahit_plot_data.rdata")
assign("megahit", fasta_frame)
load("Trans_Abyss/Trans_Abyss_plot_data.rdata")
assign("trans_abyss", fasta_frame)
load("PLASS/PLASS_plot_data.rdata")
assign("plass", fasta_frame)
load("RNASpades/RNASpades_plot_data.rdata")
assign("rna_spades", fasta_frame)
load("SOAP_DeNovo_Trans/SOAP_DeNovo_Trans_plot_data.rdata")
assign("soap_denovo_trans", fasta_frame)
load("Trinity/Trinity_plot_data.rdata")
assign("trinity", fasta_frame)

megahit %>%
  mutate(type = "megahit") -> megahit
trans_abyss %>%
  mutate(type = "trans_abyss") -> trans_abyss
plass %>%
  mutate(type = "plass") -> plass
rna_spades %>%
  mutate(type = "rna_spades") -> rna_spades
soap_denovo_trans %>%
  mutate(type = "soap_denovo_trans") -> soap_denovo_trans
trinity  %>%
  mutate(type = "trinity") -> trinity

rbind(megahit, trans_abyss, plass,
      rna_spades, soap_denovo_trans, trinity) -> contig_bins
ggplot(contig_bins, aes(x=value, fill=type)) +
  geom_histogram(color="#e9ecef", alpha=0.8, position = 'identity') +
  scale_x_log10() +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(vars(type), ncol = 1, strip.position = "right") +
  ggtitle("Assembler contig distribution, 1M reads",
          subtitle = "x-axis log10") +
  xlab("contig length") +
  scale_fill_manual(values=c("#69b3a2", "#404080", "#938BA1",
                             "#F5D491", "#B6244F", "#FF8552")) +
  labs(fill="Assembler") -> logx

ggplot(contig_bins, aes(x=value, fill=type)) +
  geom_histogram(color="#e9ecef", alpha=0.8, position = 'identity') +
  scale_y_log10() +
  scale_x_log10() +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(vars(type), ncol = 1, strip.position = "right") +
  ggtitle("Assembler contig distribution, 1M reads",
          subtitle = "x-axis + y-axis log10") +
  xlab("contig length") +
  scale_fill_manual(values=c("#69b3a2", "#404080", "#938BA1",
                             "#F5D491", "#B6244F", "#FF8552")) +
  labs(fill="Assembler") -> logboth

logx + logboth

# Plot duration -----------------------------------------------------------

data.frame(assembler = c("Megahit", "PLASS", "RNASpades",
                         "SOAP_DeNovo_Trans", "Trans_Abyss", "Trinity"),
           duration  = c(56, 187, 77, 47, 368, 3257)) -> time_data

ggplot(data = time_data, aes(x = assembler, y = duration / (60))) +
  geom_col(aes(fill = assembler)) +
  labs(title = "Assembler run duration",
       y = "Duration (Minutes)",
       x = NULL) +
  scale_fill_manual(values=c("#69b3a2", "#404080", "#938BA1",
                             "#F5D491", "#B6244F", "#FF8552")) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") -> duration
(n50 / duration) | (bowtie / contigs)
