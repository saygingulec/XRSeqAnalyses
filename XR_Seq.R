#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
sample = args[1]

setwd('results')
library(ggplot2)
theme_set(theme_minimal())

# Read Length Distribution

rl_file = paste(sample, "_read_length_distribution.txt", sep = "")

df = as.data.frame(read.delim(rl_file))
colnames(df) = c("Length", "Count")
df$Percent = prop.table(df$Count) * 100

plot = ggplot(df, aes(x = Length, y = Percent)) +
  geom_bar(stat = "identity", fill = "deepskyblue3") +
  labs(
    y = "Frequency (%)",
    x = "Lenght (nt)",
    title = sample) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  coord_cartesian(clip = "off") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste(sample, '_rld_graph.png', sep = ""), width = 5.54, height = 4.39)


# Monomer Analysis

mon_file = paste(sample, "_monomer_R_df.txt", sep = "")

df <- read.table(mon_file, header = FALSE)
colnames(df) <- c("Length", "Position", "Base", "Percentage")
df$Base <- factor(df$Base, levels = c("G", "A", "C", "T"), ordered = TRUE)
plot = ggplot(df, aes(x=Position, y=Percentage)) +
  geom_bar(stat = "identity", aes(fill = Base)) +
  labs(title = sample) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c("G" = "purple4",
                               "C" = "dodgerblue4",
                               "A" = "green4",
                               "T" = "orange")) +
  facet_grid(Length ~ .) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Nucleotide frequency (% of Total)") +
  theme(axis.text = element_text(size = 3),
        axis.title.x = (element_text(size = 10)),
        axis.title.y = (element_text(size = 10))) +
  theme(legend.title=element_text(size=10),
        legend.text=element_text(size=10),
        strip.text = element_text(size=5))

ggsave(paste(sample, '_monomer_graph.png', sep = ""), width = 5.54, height = 4.39)


# TCR Graph

nts_pat = paste(sample, ".*_NTS_rpkm_rdg_averages.txt", sep = "")
ts_pat = paste(sample, ".*_TS_rpkm_rdg_averages.txt", sep = "")
NTS_averages = list.files(pattern = nts_pat)
TS_averages = list.files(pattern = ts_pat)

files = cbind(NTS_averages, TS_averages)

for (i in 1:length(NTS_averages)) {

  NTS_name = sub("_rpkm.*", "", files[[i,1]])
  TS_name = sub("_rpkm.*", "", files[[i,2]])

  NTS_avg = read.delim(files[i, 1])$Average.Read.Count
  TS_avg = read.delim(files[i, 2])$Average.Read.Count

  df = data.frame(NTS_avg, TS_avg)
  plot = ggplot(df, aes(x = as.numeric(row.names(df)))) +
    coord_cartesian(clip = "off") +
    geom_line(aes(y = NTS_avg, color = "darkorchid1"), size = 1) +
    geom_line(aes(y = TS_avg, color = "deepskyblue3"), size = 1) +
    labs(
      y = "RPKM",
      x = "",
      title = sample) +
    theme(
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")) +
    theme(
      axis.text = element_text(size = 10),
      axis.title.y = element_text(size = 12),
      axis.title.x = element_blank(),
      plot.title = element_text(margin=margin(0,0,20,0)),
      axis.ticks = element_line(),
      plot.margin = margin(r = 10)) +
    theme(
      legend.position="bottom",
      legend.direction = "vertical",
      legend.title = element_blank()) +
    scale_color_manual(labels = c(NTS_name, TS_name),
                       values = c("darkorchid1", "deepskyblue3")) +
    scale_x_continuous(
      breaks = c(1, 25, 75, 125, 150),
      labels = c("-2 kb", "TSS", "%50", "TES", "2 kb"),
      expand = c(0, 0)) +
    scale_y_continuous(
      limits = c(0, max(max(df$TS_avg), max(df$NTS_avg)) + 1),
      expand = c(0, 0)) +
    theme(plot.title = element_text(hjust = 0.5, size = 15)) +
    geom_vline(
      xintercept = c(25, 125),
      linetype="dashed",
      color = "black", size=0.5)

  ggsave(paste(sample, '_tcr_graph.png', sep = ""), width = 5.54, height = 4.39)

}


# log2 Graph

log2_pat = paste(sample, ".*_TS_NTS_rpkm.txt", sep = "")
log2s = list.files(pattern = log2_pat)

for (i in log2s) {

  df = as.data.frame(read.delim(i))

  plot = ggplot(df, aes(x = log2)) +
    coord_cartesian(clip = "off") +
    geom_histogram(binwidth = 0.25, fill = 'gray50') +
    labs(
      y = "Frequency",
      x = "log2(TS/NTS)",
      title = sample) +
    scale_x_continuous(expand = c(0, 0),
                       breaks = seq(as.integer(min(df$log2)), as.integer(max(df$log2)), by = 2)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme(plot.title = element_text(hjust = 0.5, margin=margin(b = 10)),
          panel.grid.minor = element_blank()) +
    geom_vline(
      xintercept = c(0),
      color = "gray30", size=0.5)

  ggsave(paste(sample, '_log2_graph.png', sep = ""), width = 5.54, height = 4.39)

}
