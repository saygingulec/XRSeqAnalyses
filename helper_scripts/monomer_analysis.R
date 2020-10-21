# You only need to specify the table's name and choose the name of your output file.

library(ggplot2)

# Table's name:
df <- read.table("CEWTL1_CPD_1hR2_CCGTCC_S5_L007_R1_001_monomer_R_df.txt")

colnames(df) <- c("Length", "Position", "Base", "Percentage")
df$Base <- factor(df$Base, levels = c("G", "A", "C", "T"), ordered = TRUE)
levels(as.factor(df$Base))

plot <- ggplot(df, aes(x=Position, y=Percentage)) + 
  geom_bar(stat = "identity", aes(fill = Base)) + 
  scale_fill_manual(values = c("G" = "purple4", "C" = "dodgerblue4", "A" = "green4", "T" = "orange")) + 
  facet_grid(Length ~ .) + 
  theme(plot.title = element_text(hjust = 0.5)) + ylab("Nucleotide frequency (% of Total)") +
  theme(axis.text = element_text(size = 3), axis.title.x = (element_text(size = 10)), axis.title.y = (element_text(size = 10))) +
  theme(legend.title=element_text(size=10), legend.text=element_text(size=10), strip.text = element_text(size=5))

# Output name:
ggsave('C. Elegans WT L1 1h CPD R2.pdf', plot)

