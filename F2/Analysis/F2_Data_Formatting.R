### Formatting Data

# load files and libraries
library(DESeq2)
F2.read.count <- read.table("F2_kallisto_gene_counts.tsv", header = T, check.names = F)
rownames(F2.read.count) <- F2.read.count[,1]
F2.read.count <- F2.read.count[,-1]

# Save New Count file
save(F2.read.count, file = "F2.read.count.Rdata")

### Set up sample description

# filter based on read count, assign group, normalize, design matrix, calculate dispersion   
# set up group 

F2.read.count <- F2.read.count[,colSums(F2.read.count) > 1000000]  
dim(F2.read.count) #101040 120
Batches <- read.table("Batch.tsv", header = T, check.names = F)

F2.read.count.sample <- data.frame(colnames(F2.read.count))
colnames(F2.read.count.sample) <- "Sample"
F2.read.count.sample <- cbind(F2.read.count.sample, Batches)

save(F2.read.count.sample, file = "F2.read.count.sample.Rdata")
# filter based on read count 
F2.read.count.small <- F2.read.count[rowSums(F2.read.count > 10) >= 3,]
dim(F2.read.count.small) #63819 120
# save
save(F2.read.count.small, file = "F2.read.count.small.Rdata")

###voom transformation

dds.F2 <- DESeqDataSetFromMatrix(countData = round(F2.read.count.small), colData = F2.read.count.sample, design = ~ Batch)

vsd.F2 <- varianceStabilizingTransformation(dds.F2)
vstMat.F2<- assay(vsd.F2)
colnames(vstMat.F2) <- colnames(F2.read.count)
# save
save(vstMat.F2, file = "vstMat.F2.Rdata")
