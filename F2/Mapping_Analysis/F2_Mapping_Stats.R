### Mapping stats anaylisis
### setwd("~/Documents/KIAT/F2")

### Load in stats and libraries
library(ggplot2)
library(reshape2)
library(scales)
star_stats <- read.csv("modified.modified.F2.stats.csv")
star_stats$Batch <- sub("1", "a", star_stats$Batch)
star_stats$Batch <- sub("2", "b", star_stats$Batch)

### Summary Barplot of mapping results

s_star_stats <- star_stats[,c(1,8,10)]
s_star_stats$Unmapped <- 100 -s_star_stats$Percent_Unique_Mapped - s_star_stats$Percent_Multi_Mapped
colnames(s_star_stats) <- c("Sample_ID","Uniquely Mapped", "Multi Mapped", "Unmapped")
m_star_stats <- melt(s_star_stats)
m_star_stats <- m_star_stats[order(m_star_stats$variable, m_star_stats$value),]
m_star_stats$variable <- factor(m_star_stats$variable, levels = rev(levels(m_star_stats$variable)))
m_star_stats$Sample_ID <- factor(m_star_stats$Sample_ID, levels = m_star_stats$Sample_ID)

ggplot(m_star_stats, aes(factor(Sample_ID), value, fill = variable, order = as.numeric(variable))) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x=element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  scale_y_continuous(breaks = seq(0, 100, by = 5)) +
  ggtitle("Overall Mapping Results") + xlab("Sample ID") +
  ylab("Percentage") + labs(fill = "Mapping")

#ggsave("Plots/Overall_Mapping_Results.png", width = 17, height = 10)

### Barplots total reads
#### Total Reads

ggplot(star_stats, aes(reorder(Sample_ID, TotalReads), TotalReads)) + 
  geom_bar(stat = "identity", fill = 'blue') +
  ggtitle("Barplot of Total Reads") +
  xlab("Sample") +
  ylab("Number of Reads") +
  theme(axis.text.x=element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5))

#ggsave("Plots/Total_Reads_Barplot.png", width = 17)
  
### Total Reads + Type

s_star_stats <- star_stats[,c(1,3)]
s_star_stats$LQ.reads <- star_stats$TotalReads - star_stats$HQ.reads
colnames(s_star_stats) <- c('Sample_ID', 'High Quality', 'Low Quality')
m_star_stats <- melt(s_star_stats)

ggplot(m_star_stats, aes(reorder(Sample_ID, value), value, fill = variable)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x=element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  ggtitle("Barplot of Total Reads and Type") + xlab("Sample ID") +
  ylab("Number of Reads") + labs(fill = "Quality")

#ggsave("Plots/Total_Reads_Plus_Type_Barplot.png", width = 17)


### Mapped, unmapped, trimmed off

s_star_stats <- data.frame(star_stats[,1])
s_star_stats$Mapped <- star_stats[,7] + star_stats[,9]
s_star_stats$Unmapped <- star_stats[,3] - s_star_stats[,2]
s_star_stats$Trimmed_Off <- star_stats[,2] - star_stats[,3]
colnames(s_star_stats) <- c('Sample_ID', 'Mapped', 'Unmapped', 'Trimmed Off')
m_star_stats <- melt(s_star_stats)

ggplot(m_star_stats, aes(reorder(Sample_ID, value), value, fill = variable)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x=element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  ggtitle("Barplot of Mapping Results") + xlab("Sample ID") +
  ylab("Number of Reads") + labs(fill = "Type")

#ggsave("Plots/Mapping_Barplot.png", width = 17)

### Unique, Multi, Unmapped, Trimmed Off

s_star_stats <- data.frame(star_stats[,1])
s_star_stats$Unique <- star_stats[,7] 
s_star_stats$Multi <- star_stats[,9]
s_star_stats$Unmapped <- star_stats[,3] - s_star_stats[,2] - s_star_stats[,3]
s_star_stats$Trimmed_Off <- star_stats[,2] - star_stats[,3]
colnames(s_star_stats) <- c('Sample_ID', 'Unique', 'Multi', 'Unmapped', 'Trimmed Off')
m_star_stats <- melt(s_star_stats)

ggplot(m_star_stats, aes(reorder(Sample_ID, value), value, fill = variable)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x=element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  ggtitle("Barplot of Mapping Results") + xlab("Sample ID") +
  ylab("Number of Reads") + labs(fill = "Type")

#ggsave("Plots/Unique_Mapping_Barplot.png", width = 17)

### High Quality Reads

ggplot(star_stats, aes(reorder(Sample_ID, HQ.reads), HQ.reads)) + 
  geom_bar(stat = "identity", fill = "red") +
  ggtitle("Barplot of Total HQ Reads") +
  xlab("Sample") +
  ylab("Number of Reads") +
  theme(axis.text.x=element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5))

#ggsave("Plots/High_Quality_Reads_Barplot.png", width = 17)

### Average Total High Quality Reads

ggplot(star_stats, aes(HQ.reads)) + 
  geom_histogram(bins = 40, color = "blue", fill = "blue") +
  ggtitle("Histogram of High Quality Reads") +
  xlab("Number of high quality reads") +
  ylab("Count") +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 1.8e7, y = 21, label = paste("Average number of reads", round(mean(star_stats$HQ.reads))))

#ggsave("Plots/Histogram_HQ_Reads.png", width = 10)

### Average Mapped Reads

s_star_stats <- data.frame(star_stats[,1])
s_star_stats$Mapped <- star_stats[,7] + star_stats[,9]

ggplot(s_star_stats, aes(Mapped)) + 
  geom_histogram(bins = 40, color = "blue", fill = "blue") +
  ggtitle("Histogram of Mapped Reads") +
  xlab("Number of mapped reads") +
  ylab("Count") +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 1.3e7, y = 13, label = paste("Average number of mapped reads", round(mean(s_star_stats$Mapped))))

#ggsave("Plots/Histogram_Mapped_Reads.png", width = 10)

## Average of mapping percentage results
Unique <- round(sum(star_stats$Percent_Unique_Mapped)/nrow(star_stats),2)
Multi <- round(sum(star_stats$Percent_Multi_Mapped)/nrow(star_stats),2)
Unmapped <- 100 - Unique - Multi

total <- melt(data.frame(Unique,Multi,Unmapped))

ggplot(total, aes(reorder(variable, value), value, fill = variable)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x=element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  scale_y_continuous(breaks = seq(0, 100, by = 5)) +
  ggtitle("Average of Mapping Results") + xlab("Mapping") +
  ylab("Percentage") + labs(fill = "Mapping")

#ggsave("Plots/Average_Overall_Mapping_Results.png")

### Histogram of Percent Unique

ggplot(star_stats, aes(Percent_Unique_Mapped)) + 
  geom_histogram(bins = 20, color = "blue", fill = "blue") +
  ggtitle("Histogram of Uniquely Mapped Reads") +
  xlab("Percent of uniquely mapped reads") +
  ylab("Count") +
  theme(plot.title = element_text(hjust = 0.5))

#ggsave("Plots/Percent_Unique_Mapping_Histogram.png", width = 7)

### Total reads by Batch

ggplot(star_stats, aes(reorder(Sample_ID, TotalReads), TotalReads, fill = Batch)) + 
  geom_bar(stat = "identity") +
  ggtitle("Barplot of Total Reads") +
  xlab("Sample") +
  ylab("Number of Reads") +
  theme(axis.text.x=element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5))


### Histogram of Percent Unique by batch

ggplot(star_stats, aes(Percent_Unique_Mapped, fill = Batch)) + 
  geom_histogram(bins = 20) +
  ggtitle("Histogram of Uniquely Mapped Reads") +
  xlab("Percent of uniquely mapped reads") +
  ylab("Count") +
  theme(plot.title = element_text(hjust = 0.5))

#ggsave("Batch_Plots/Percent_Unique_Histogram.png", width = 6)

### Histogram of All multimapped reads By Batch

s_star_stats <- data.frame(star_stats$Sample_ID)
s_star_stats$Total_Percent_Multi_Mapped <- star_stats$Percent_Multi_Mapped + star_stats$Percent_Too_Many_Multi_Mapped
s_star_stats$Batch <- star_stats$Batch

ggplot(s_star_stats, aes(Total_Percent_Multi_Mapped, fill = Batch)) + 
  geom_histogram(bins = 20) +
  ggtitle("Histogram of All Multi-mapped Reads") +
  xlab("Percent of Multi-mapped Reads") +
  ylab("Count") +
  theme(plot.title = element_text(hjust = 0.5))

#ggsave("Batch_Plots/All_Multi-Mapped_Histogram.png", width = 7)

### Histogram of multimapped reads (<10 loci matches) By Batch

ggplot(star_stats, aes(Percent_Multi_Mapped, fill = Batch)) + 
  geom_histogram(bins = 20) +
  ggtitle("Histogram of Multi-mapped Reads (<10 loci)") +
  xlab("Percent of multi-mapped reads") +
  ylab("Count") +
  theme(plot.title = element_text(hjust = 0.5))

#ggsave("Batch_Plots/Multi-Mapped_Histogram.png", width = 7)

### Histogram of too many multi (>= 10 loci matches) By Batch

ggplot(star_stats, aes(Percent_Too_Many_Multi_Mapped, fill = Batch)) + 
  geom_histogram(bins = 20) +
  ggtitle("Histogram of Mapped to too many Reads (>=10 loci") +
  xlab("Percent of mapped to too many reads") +
  ylab("Count") +
  theme(plot.title = element_text(hjust = 0.5))

#ggsave("Batch_Plots/Too_Many_Multi_Histogram.png", width = 7)

### Histogram of True Unmapped By Batch

s_star_stats <- data.frame(star_stats$Sample_ID)
s_star_stats$Unmapped <- star_stats$Percent_Unmapped_Mismatches + star_stats$Percent_Unmapped_Too_Short + star_stats$Percent_Unmapped_Other
s_star_stats$Batch <- star_stats$Batch

ggplot(s_star_stats, aes(Unmapped, fill = Batch)) + 
  geom_histogram(bins = 20) +
  ggtitle("Histogram of Unmapped Reads") +
  xlab("Percent of unmapped reads") +
  ylab("Count") +
  theme(plot.title = element_text(hjust = 0.5))

#ggsave("Batch_Plots/Unmapped_Histogram.png", width = 7)

### Percent High Quality Reads By Batch

ggplot(star_stats, aes(HQ_percent, fill = Batch)) + 
  geom_histogram() +
  ggtitle("Histogram of High Quality Reads") +
  xlab("Percent of high quality reads") +
  ylab("Count") +
  theme(plot.title = element_text(hjust = 0.5))

#ggsave("Batch_Plots/High_Quality_Reads_Histogram.png", width= 7)

### High Quality Reads by Batch

ggplot(star_stats, aes(reorder(Sample_ID, HQ.reads), HQ.reads, fill = Batch)) + 
  geom_bar(stat = "identity") +
  ggtitle("Barplot of Total HQ Reads") +
  xlab("Sample") +
  ylab("Number of Reads") +
  theme(axis.text.x=element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5))

## Average of mapping percentage results by batch
a.star_stats <- subset(star_stats, star_stats$Batch == "a")
b.star_stats <- subset(star_stats, star_stats$Batch == "b")
Unique <- round(sum(a.star_stats$Percent_Unique_Mapped)/nrow(a.star_stats),2)
Multi <- round(sum(a.star_stats$Percent_Multi_Mapped)/nrow(a.star_stats),2)
Unmapped <- 100 - Unique - Multi
Batch <- "a"
a.total <- melt(data.frame(Unique,Multi,Unmapped,Batch))
Unique <- round(sum(b.star_stats$Percent_Unique_Mapped)/nrow(b.star_stats),2)
Multi <- round(sum(b.star_stats$Percent_Multi_Mapped)/nrow(b.star_stats),2)
Unmapped <- 100 - Unique - Multi
Batch <- "b"
b.total <- melt(data.frame(Unique,Multi,Unmapped,Batch))
total <- rbind(a.total,b.total)
positions <- c("Unmapped","Multi","Unique")

ggplot(total, aes(variable, value, fill = variable)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x=element_text(angle=90, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  ggtitle("Average of Mapping Results") + xlab("Mapping") +
  ylab("Percentage") + labs(fill = "Mapping") +
  facet_grid(~Batch) +
  scale_x_discrete(limits = positions)

#ggsave("Batch_Plots/Average_Mapping_Rates_Batch.png")
