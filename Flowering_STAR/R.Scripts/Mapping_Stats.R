### looking at mapping stats
library(reshape2)
library(ggplot2)
stats <- read.delim2("mapping.stats", header = TRUE, sep = "\t", col.names = c("Unique","Multi","Unmapped",""))
stats <- stats[,1:3]

###Plotting
melt.stats <- melt(t(stats))
melt.stats$value <- as.numeric(as.character(melt.stats$value))
colnames(melt.stats) <- c("Mapping","Sample","Percentage")
# plot and facet by Condition
ggplot(data=melt.stats, aes(Mapping,Percentage, fill=Mapping)) + facet_wrap(~ Sample, scales = "free_y") + geom_bar(stat = "identity")
ggsave("Mapping.barplot.png", width = 20, height = 8)
