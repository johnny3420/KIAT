###Giving descriptions to DE genes
###Load data
anno <- read.delim("Brassica_napus_IPR.withdescription.gz", header = F, sep = "\t")
DEgene.gt.100bp.paired <- read.csv("Flowering.DE_genes/DEgene.gt.100bp.paired.csv")
DEgene.gt.50bp.paired <- read.csv("Flowering.DE_genes/DEgene.gt.50bp.paired.csv")
DEgene.gt.100bp.single <- read.csv("Flowering.DE_genes/DEgene.gt.100bp.single.csv")
DEgene.gt.50bp.single <- read.csv("Flowering.DE_genes/DEgene.gt.50bp.single.csv")

###Format data
anno <- anno[,c(1,3)]
p100 <- as.data.frame(sort(DEgene.gt.100bp.paired[,1]))
ap100 <- merge(p100,anno, by.x = 1, by.y = 1, all.x = T)
p50 <- as.data.frame(sort(DEgene.gt.50bp.paired[,1]))
ap50 <- merge(p50,anno, by.x = 1, by.y = 1, all.x = T)
s100 <- as.data.frame(sort(DEgene.gt.100bp.single[,1]))
as100 <- merge(s100,anno, by.x = 1, by.y = 1, all.x = T)
s50 <- as.data.frame(sort(DEgene.gt.50bp.single[,1]))
as50 <- merge(s50,anno, by.x = 1, by.y = 1, all.x = T)

###Write new tables
write.table(ap100, "Flowering.DE_genes/100bp_paired_end_DE_genes_and_descriptions.tab",sep = "\t", col.names = F, row.names = F)
write.table(ap50, "Flowering.DE_genes/50bp_paired_end_DE_genes_and_descriptions.tab",sep = "\t", col.names = F, row.names = F)
write.table(as100, "Flowering.DE_genes/100bp_single_end_DE_genes_and_descriptions.tab",sep = "\t", col.names = F, row.names = F)
write.table(as50, "Flowering.DE_genes/50bp_single_end_DE_genes_and_descriptions.tab",sep = "\t", col.names = F, row.names = F)

