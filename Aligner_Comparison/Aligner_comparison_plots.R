###Aligner comparision, quick visualization
library(ggplot2)
library(reshape2)
stats <- read.delim("Aligner_Comparison/Aligner_comparison.tab", header = T, sep = "\t")
melt.stats <- melt(stats)
colnames(melt.stats) <- c("Alignment","Aligner","Mapped")
ggplot(data=melt.stats, aes(Aligner,Mapped, fill=Aligner)) + 
  facet_wrap("Alignment") + 
  geom_bar(stat = "identity")
#ggsave("Aligner_comparison_barplot.png")

m.stats <- read.delim("Aligner_Comparison/modified.aligner.comparison.tab", header = T, sep = "\t")
melt.m.stats <- melt(m.stats)
colnames(melt.m.stats) <- c("Sample","Aligner","Mapping","Percentage")

ggplot(data=melt.m.stats, aes(Mapping,Percentage, fill=Mapping)) + 
  facet_grid(Sample ~ Aligner) + 
  geom_bar(stat = "identity")

ggsave("modified.aligner_comparison_barplot.png")
