### MDS and Clustering

# load files and Libraries
library(edgeR)
library(ggplot2)
library(DESeq2)
library(ggdendro)
load("vstMat.flowering.Rdata")
load("flowering.read.count.small.Rdata")

# Normalize
dge.data.flowering <- DGEList(counts=flowering.read.count.small, group=flowering.read.count.sample$group)
dge.data.flowering <- calcNormFactors(dge.data.flowering, method = "TMM") 
dge.data.flowering$sample
hist(dge.data.flowering$samples$norm.factors)

# MDS to check sample seperation, also using dengrogram.  
mds <- plotMDS(dge.data.flowering, method = "bcv",labels = dge.data.flowering$samples$group)

x <- as.data.frame(mds$x)
y <- as.data.frame(mds$y)
distance_matrix <- merge(x, y, by="row.names")
distance_matrix$group <- gsub("(Da-Ae|Da-Ol-1)(_)(Young|flowering|early-silique|late-silique|bolting)(_)(100bp|50bp)(_)(paired|single)(_)(1|2|3)","\\5\\7\\1",distance_matrix$Row.names)
distance_matrix$gt <- gsub("(Da-Ae|Da-Ol-1)(_)(Young|flowering|early-silique|late-silique|bolting)(_)(100bp|50bp)(_)(paired|single)(_)(1|2|3)","\\1",distance_matrix$Row.names)
distance_matrix$length <- gsub("(Da-Ae|Da-Ol-1)(_)(Young|flowering|early-silique|late-silique|bolting)(_)(100bp|50bp)(_)(paired|single)(_)(1|2|3)","\\5",distance_matrix$Row.names)
distance_matrix$pairing <- gsub("(Da-Ae|Da-Ol-1)(_)(Young|flowering|early-silique|late-silique|bolting)(_)(100bp|50bp)(_)(paired|single)(_)(1|2|3)","\\7",distance_matrix$Row.names)

colnames(distance_matrix) <- c("lib","x","y","group","gt","length","pairing")
head(distance_matrix)

# making color MDS figure 
p.mds <- ggplot(data = distance_matrix)
p.mds <- p.mds + geom_point(aes(x, y, color=factor(length)), size=3) + geom_text(aes(x,y,label=group),hjust=0, vjust=0)
p.mds <- p.mds + labs(y = "BCV distance 2", x="BCV distance 1")
# p.mds <- p.mds + facet_grid(~gt)
p.mds

# Save Plot
ggsave("MDS.png", width = 11, height = 8)

# clustering to check sample seperation 


hc.flowering <- hclust(dist(t(vstMat.flowering)))

ggdata.flowering <- dendro_data(as.dendrogram(hc.flowering))
ggdata.flowering$labels$gt <- gsub("(Da-Ae|Da-Ol-1)(_)(Young|flowering|early-silique|late-silique|bolting)(_)(100bp|50bp)(_)(paired|single)(_)(1|2|3)","\\1",ggdata.flowering$labels$label) 
ggdata.flowering$labels$pairing <- gsub("(Da-Ae|Da-Ol-1)(_)(Young|flowering|early-silique|late-silique|bolting)(_)(100bp|50bp)(_)(paired|single)(_)(1|2|3)","\\7",ggdata.flowering$labels$label)
ggdata.flowering$labels$length <- gsub("(Da-Ae|Da-Ol-1)(_)(Young|flowering|early-silique|late-silique|bolting)(_)(100bp|50bp)(_)(paired|single)(_)(1|2|3)","\\5",ggdata.flowering$labels$label)
ggdata.flowering$labels$group <- paste(ggdata.flowering$labels$length, ggdata.flowering$labels$pairing, ggdata.flowering$labels$gt, sep = "_")
ggdata.flowering$labels

# start ggplot here 
p.dengro <- ggplot(data = segment(ggdata.flowering))
p.dendro <- p.dengro + geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + theme_dendro()
p.dendro <- p.dendro + geom_text(data = label(ggdata.flowering), aes(x = x, y = y, label = label, hjust=0, color=group)) + coord_flip() + scale_y_reverse(expand=c(0.2, 0))
p.dendro

# Save Plot
ggsave("clustering_new.png", width = 20, height = 8)
