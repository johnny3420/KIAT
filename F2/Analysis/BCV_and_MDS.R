### MDS and Clustering

# load files and Libraries
library(edgeR)
library(ggplot2)
library(DESeq2)
library(ggdendro)
load("vstMat.F2.Rdata")
load("F2.read.count.small.Rdata")
load("F2.read.count.sample.Rdata")

# Normalize
dge.data.F2 <- DGEList(counts=F2.read.count.small, group = F2.read.count.sample$Batch)
dge.data.F2 <- calcNormFactors(dge.data.F2, method = "TMM") 
dge.data.F2$sample
hist(dge.data.F2$samples$norm.factors)

# MDS to check sample seperation, also using dengrogram.  
mds <- plotMDS(dge.data.F2, method = "bcv",labels = dge.data.F2$samples$group)

x <- as.data.frame(mds$x)
y <- as.data.frame(mds$y)
distance_matrix <- merge(x, y, by="row.names")
distance_matrix$group <- dge.data.F2$samples$group
colnames(distance_matrix) <- c("lib","x","y","Batch")
head(distance_matrix)

# making color MDS figure 
p.mds <- ggplot(data = distance_matrix)
p.mds <- p.mds + geom_point(aes(x, y, color=factor(Batch)), size=3)
p.mds <- p.mds + labs(y = "BCV distance 2", x="BCV distance 1")
p.mds

# Save Plot
ggsave("MDS.png", height = 10, width = 10)

# clustering to check sample seperation 


hc.F2 <- hclust(dist(t(vstMat.F2)))

ggdata.F2 <- dendro_data(as.dendrogram(hc.F2))
ggdata.F2$labels$Batch <- F2.read.count.sample$Batch
ggdata.F2$labels

# start ggplot here 
p.dengro <- ggplot(data = segment(ggdata.F2))
p.dendro <- p.dengro + geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + theme_dendro()
p.dendro <- p.dendro + geom_text(data = label(ggdata.F2), aes(x = x, y = y, label = label, hjust=0, color=Batch)) + coord_flip() + scale_y_reverse(expand=c(0.1, 0))
p.dendro

# Save Plot
ggsave("clustering.png", width = 20, height = 15)
