###Giving descriptions to DE genes
###Load data
library(plyr)
anno <- read.delim("Brassica_napus_IPR.withdescription.gz", header = F, sep = "\t")
DEgene.gt.100bp.paired <- read.csv("R_Analysis/DEgene.gt.100bp.paired.csv")
DEgene.gt.50bp.paired <- read.csv("R_Analysis/DEgene.gt.50bp.paired.csv")
DEgene.gt.100bp.single <- read.csv("R_Analysis/DEgene.gt.100bp.single.csv")
DEgene.gt.50bp.single <- read.csv("R_Analysis/DEgene.gt.50bp.single.csv")

###Format data
anno <- anno[,c(1,3)]
condensed_anno <- ddply(anno, .(V1), summarize, V2=paste(V3,collapse=", "))
colnames(condensed_anno) <- c("Gene","Description")
p100 <- as.data.frame(sort(DEgene.gt.100bp.paired[,1]))
ap100 <- merge(p100,condensed_anno, by.x = 1, by.y = 1, all.x = T)
p50 <- as.data.frame(sort(DEgene.gt.50bp.paired[,1]))
ap50 <- merge(p50,condensed_anno, by.x = 1, by.y = 1, all.x = T)
s100 <- as.data.frame(sort(DEgene.gt.100bp.single[,1]))
as100 <- merge(s100,condensed_anno, by.x = 1, by.y = 1, all.x = T)
s50 <- as.data.frame(sort(DEgene.gt.50bp.single[,1]))
as50 <- merge(s50,condensed_anno, by.x = 1, by.y = 1, all.x = T)

###Write new tables
write.table(ap100, "R_Analysis/100bp_paired_end_DE_genes_and_descriptions.tab",sep = "\t", col.names = F, row.names = F)
write.table(ap50, "R_Analysis/50bp_paired_end_DE_genes_and_descriptions.tab",sep = "\t", col.names = F, row.names = F)
write.table(as100, "R_Analysis/100bp_single_end_DE_genes_and_descriptions.tab",sep = "\t", col.names = F, row.names = F)
write.table(as50, "R_Analysis/50bp_single_end_DE_genes_and_descriptions.tab",sep = "\t", col.names = F, row.names = F)

