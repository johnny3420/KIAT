### Formatting Data

# load files and libraries
library(DESeq2)
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