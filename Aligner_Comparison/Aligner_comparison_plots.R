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
