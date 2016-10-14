### looking at mapping stats
library(reshape2)
library(ggplot2)
stats <- read.delim2("Flowering_STAR/Mapping.stats", header = TRUE, sep = "\t", row.names = NULL, col.names = c("Unique","Multi","Unmapped","Reads"))
row.names(stats) <- stats[,1]
stats <- stats[,c(-1,-5)]
### Formating data frame
melt.stats <- melt(t(stats))
melt.stats$value <- as.numeric(as.character(melt.stats$value))
colnames(melt.stats) <- c("Mapping","Sample","Percentage")
melt.stats$Size <- sub("\\..*$","", melt.stats$Sample)
melt.stats$Pairing <- sub("\\..*$","", sub(".*bp\\.","", melt.stats$Sample))
melt.stats$Gt <- sub("All1_Gae_2|All1_Gae_3","Da-Ol-1",sub(".*\\.","",melt.stats$Sample))
melt.stats$Gt <-sub("2","Da-Ol-1",sub("6|Ae_Gae_2|Ae_Gae_3","Da-Ae",melt.stats$Gt))
melt.stats$Replicate <- sub(".*\\.","",melt.stats$Sample)
# Plotting
ggplot(data=melt.stats, aes(Mapping,Percentage, fill=Mapping)) + facet_wrap(~Sample, scales = "free_y") + geom_bar(stat = "identity")
ggsave("R_Analysis/Mapping.barplot.png")
# plot and facet by Pairing and Size
ggplot(data=melt.stats, aes(Sample,Percentage, fill=Mapping)) + facet_grid(Size~Pairing, scales = "free_x") + geom_bar(stat = "identity") + theme(axis.text.x=element_text(angle=90, hjust=1))
ggsave("R_Analysis/Mapping.barplot.2.png")
# plot and facet by Genotype
ggplot(data=melt.stats, aes(Sample,Percentage, fill=Mapping)) + facet_grid(~Gt, scales = "free_x") + geom_bar(stat = "identity") + theme(axis.text.x=element_text(angle=90, hjust=1))
ggsave("R_Analysis/Mapping.barplot.3.png")
# plot and facet by Sample
ggplot(data=melt.stats, aes(Mapping,Percentage, fill=Mapping)) + facet_grid(~Sample, scales = "free_x") + 
  geom_bar(stat = "identity") + theme(axis.text.x=element_text(angle=90, hjust=1)) + ggtitle("Mapping Rates")
# plot and facet by replicate
ggplot(data=melt.stats, aes(Sample,Percentage, fill=Mapping)) + facet_grid(~Replicate, scales = "free_x") + 
  geom_bar(stat = "identity") + theme(axis.text.x=element_text(angle=90, hjust=1)) + ggtitle("Mapping Rates")
ggsave("R_Analysis/Mapping.barplot.4.png")
##
melt.stats$Condition <- sub("\\.All1_Gae_2|\\.All1_Gae_3|\\.2|\\.6|\\.Ae_Gae_2|\\.Ae_Gae_3","",melt.stats$Sample)
ggplot(data=melt.stats, aes(Mapping,Percentage/6, fill=Mapping)) + ylab("Percent Mapped") + facet_grid(~Condition, scales = "free_x") + 
  geom_bar(stat = "identity") + theme(axis.text.x=element_text(angle=90, hjust=1)) + ggtitle("Overall Mapping Rates")
ggsave("R_Analysis/Mapping.barplot.5.png", width = 8)
