# load necessary libs & functions 
library(edgeR)
library(ggplot2)
library(DESeq2)
library(ggdendro)
library(reshape2)
### Formatting Data
flowering.read.count <- read.table("Counts/flowering.star.read.count.tsv", header = T, check.names = F)
rownames(flowering.read.count) <- flowering.read.count[,1]
flowering.read.count <- flowering.read.count[,-1]

# Sample description
sample_des <- read.csv("parent_summary_corrected.csv")
sorted_sample_des <- sample_des[order(sample_des$Stage),]
flowering_sample_des <-sorted_sample_des[10:15,]
flowering_sample_des <- flowering_sample_des[order(flowering_sample_des$SampleID),]
flowering_sample_des <- flowering_sample_des[rep(seq_len(nrow(flowering_sample_des)), 4),]
flowering_sample_des$Length <- c(rep("50bp",12), rep("100bp",12))
flowering_sample_des$Pairing <- c(rep("paired",6), rep("single",6))
flowering_sample_des <- flowering_sample_des[,c(4:7, 23:24)]

flowering_sample_ID <- paste(flowering_sample_des$Cultivar, flowering_sample_des$Stage, flowering_sample_des$Length, flowering_sample_des$Pairing, flowering_sample_des$rep, sep = "_")
flowering_sample_ID

# Replace sample ID
colnames(flowering.read.count) <- flowering_sample_ID
head(flowering.read.count)

# Save New Count file
save(flowering.read.count, file = "flowering.read.count.Rdata")

### Set up sample description


# filter based on read count, assign group, normalize, design matrix, calculate dispersion   
# set up group 

load("flowering.read.count.Rdata")
flowering.read.count <- flowering.read.count[,colSums(flowering.read.count) > 1000000]  
dim(flowering.read.count) #101040 24
flowering.read.count.sample<-data.frame(file=colnames(flowering.read.count),
                                         batch=factor(gsub("(Da-Ae|Da-Ol-1)(_)(Young|flowering|early-silique|late-silique|bolting)(_)(100bp|50bp)(_)(paired|single)(_)(1|2|3)","\\9",colnames(flowering.read.count))),  
                                         genotype=factor(gsub("(Da-Ae|Da-Ol-1)(_)(Young|flowering|early-silique|late-silique|bolting)(_)(100bp|50bp)(_)(paired|single)(_)(1|2|3)","\\1",colnames(flowering.read.count))),
                                         length=factor(gsub("(Da-Ae|Da-Ol-1)(_)(Young|flowering|early-silique|late-silique|bolting)(_)(100bp|50bp)(_)(paired|single)(_)(1|2|3)", "\\5",colnames(flowering.read.count))),
                                         pairing=factor(gsub("(Da-Ae|Da-Ol-1)(_)(Young|flowering|early-silique|late-silique|bolting)(_)(100bp|50bp)(_)(paired|single)(_)(1|2|3)", "\\7",colnames(flowering.read.count))),
                                         stage=factor(gsub("(Da-Ae|Da-Ol-1)(_)(Young|flowering|early-silique|late-silique|bolting)(_)(100bp|50bp)(_)(paired|single)(_)(1|2|3)","\\3",colnames(flowering.read.count))),  
                                         group=factor(gsub("(Da-Ae|Da-Ol-1)(_)(Young|flowering|early-silique|late-silique|bolting)(_)(100bp|50bp)(_)(paired|single)(_)(1|2|3)","\\5\\7\\1",colnames(flowering.read.count)))
)

flowering.read.count.sample
#ftable(flowering.read.count.sample,row.vars="length",col.vars=c("batch","genotype"))

# filter based on read count 
flowering.read.count.small <- flowering.read.count[rowSums(flowering.read.count > 10) >= 3,]
dim(flowering.read.count.small) #56498 24
# save
save(flowering.read.count.small, file = "flowering.read.count.small.Rdata")

###voom transformation

load("flowering.read.count.small.Rdata")

dds.flowering <- DESeqDataSetFromMatrix(countData = round(flowering.read.count.small), colData = flowering.read.count.sample, design = ~ batch + genotype)

vsd.flowering <- varianceStabilizingTransformation(dds.flowering)
vstMat.flowering <- assay(vsd.flowering)
colnames(vstMat.flowering) <- colnames(flowering.read.count)
# save
save(vstMat.flowering, file = "vstMat.flowering.Rdata")

### MDS, clustering, and expression analysis

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
ggsave("clustering_new.png", width = 15, height = 8)

# pairwise comparison using GLM model 
# design matrix
design.flowering <- model.matrix(~0+group, data = flowering.read.count.sample)
dge.data.flowering <- estimateGLMCommonDisp(dge.data.flowering, design.flowering,verbose = TRUE)
dge.data.flowering <- estimateGLMTrendedDisp(dge.data.flowering,design.flowering)
dge.data.flowering <- estimateGLMTagwiseDisp(dge.data.flowering,design.flowering)
plotBCV(dge.data.flowering)
