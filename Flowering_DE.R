### Comparing different conditions

#Load files and library
library(edgeR)
library(reshape2)
library(ggplot2)
load("flowering.read.count.small.Rdata")
load("flowering.read.count.sample.Rdata")

# Create DGElist and Normalize
dge.data.flowering <- DGEList(counts=flowering.read.count.small, group=flowering.read.count.sample$group)
dge.data.flowering <- calcNormFactors(dge.data.flowering, method = "TMM") 

# Setting up design
design.flowering <- model.matrix(~0+group, data=dge.data.flowering$samples)
colnames(design.flowering) <- levels(dge.data.flowering$samples$group)
design.flowering

# estimating dispersions and create fit
dge.data.flowering <- estimateGLMCommonDisp(dge.data.flowering, design.flowering,verbose = TRUE) # Disp = 0.1588 , BCV = 0.3985
dge.data.flowering <- estimateGLMTrendedDisp(dge.data.flowering,design.flowering)
dge.data.flowering <- estimateGLMTagwiseDisp(dge.data.flowering,design.flowering)
design.flowering
fit <- glmFit(dge.data.flowering, design.flowering)

## Finding differentially expressed genes within genotype and pairing based on read length
# Da-Ae paired
DA.P.100vs50.lrt <-glmLRT(fit, contrast = c(1,0,0,0,-1,0,0,0))
DA.P.100vs50 <- summary(decideTestsDGE(DA.P.100vs50.lrt,p.value=0.01))
# Da-Ol-1 paired
DO.P.100vs50.lrt <-glmLRT(fit, contrast = c(0,1,0,0,0,-1,0,0))
DO.P.100vs50 <- summary(decideTestsDGE(DO.P.100vs50.lrt,p.value=0.01))
# Da-Ae single
DA.S.100vs50.lrt <-glmLRT(fit, contrast = c(0,0,1,0,0,0,-1,0))
DA.S.100vs50 <- summary(decideTestsDGE(DA.S.100vs50.lrt,p.value=0.01))
# Da-Ol-1 single
DO.S.100vs50.lrt <-glmLRT(fit, contrast = c(0,0,0,1,0,0,0,-1))
DO.S.100vs50 <- summary(decideTestsDGE(DO.S.100vs50.lrt,p.value=0.01))

length <- cbind(DA.P.100vs50,DO.P.100vs50,DA.S.100vs50,DO.S.100vs50)
colnames(length) <- c("DA.P.100vs50","DO.P.100vs50","DA.S.100vs50","DO.S.100vs50")

## Finding differentially expressed genes within genotype and read length, based on pairing

# Da-Ae 100bp
DA.100.PvsS.lrt <-glmLRT(fit, contrast = c(1,0,-1,0,0,0,0,0))
DA.100.PvsS <- summary(decideTestsDGE(DA.100.PvsS.lrt,p.value=0.01))
# Da-Ol-1 100bp
DO.100.PvsS.lrt <-glmLRT(fit, contrast = c(0,1,0,-1,0,0,0,0))
DO.100.PvsS <- summary(decideTestsDGE(DO.100.PvsS.lrt,p.value=0.01))
# Da-Ae 50bp
DA.50.PvsS.lrt <-glmLRT(fit, contrast = c(0,0,0,0,1,0,-1,0))
DA.50.PvsS <- summary(decideTestsDGE(DA.50.PvsS.lrt,p.value=0.01))
# Da-Ol-1 50bp
DO.50.PvsS.lrt <-glmLRT(fit, contrast = c(0,0,0,0,0,1,0,-1))
DO.50.PvsS <- summary(decideTestsDGE(DO.50.PvsS.lrt,p.value=0.01))

pairing <- cbind(DA.100.PvsS,DO.100.PvsS,DA.50.PvsS,DO.50.PvsS)
colnames(pairing) <- c("DA.100.PvsS","DO.100.PvsS","DA.50.PvsS","DO.50.PvsS")

## Finding differentially expressed genes between genotypes, constant pairing and length

# 100bp paired
DAvsDO.100bp.paired.lrt <-glmLRT(fit, contrast = c(1,-1,0,0,0,0,0,0))
DAvsDO.100bp.paired <- summary(decideTestsDGE(DAvsDO.100bp.paired.lrt,p.value=0.01))
# 100bp single
DAvsDO.100bp.single.lrt <-glmLRT(fit, contrast = c(0,0,1,-1,0,0,0,0))
DAvsDO.100bp.single <- summary(decideTestsDGE(DAvsDO.100bp.single.lrt,p.value=0.01))
# 50bp paired
DAvsDO.50bp.paired.lrt <-glmLRT(fit, contrast = c(0,0,0,0,1,-1,0,0))
DAvsDO.50bp.paired <- summary(decideTestsDGE(DAvsDO.50bp.paired.lrt,p.value=0.01))
# 50bp single
DAvsDO.50bp.single.lrt <-glmLRT(fit, contrast = c(0,0,0,0,0,0,1,-1))
DAvsDO.50bp.single <- summary(decideTestsDGE(DAvsDO.50bp.single.lrt,p.value=0.01))

genotype <- cbind(DAvsDO.100bp.paired,DAvsDO.100bp.single,DAvsDO.50bp.paired,DAvsDO.50bp.single)
colnames(genotype) <- c('DAvsDO.100bp.paired','DAvsDO.100bp.single','DAvsDO.50bp.paired','DAvsDO.50bp.single')

Flowering.DE <- as.data.frame(cbind(genotype,length,pairing))
Flowering.DE[4,] <- colSums(Flowering.DE)
for (i in 1:ncol(Flowering.DE)){
  Flowering.DE[5,i] <- (Flowering.DE[1,i] + Flowering.DE[3,i])
  Flowering.DE[6,i] <- sprintf("%.2f", round(((Flowering.DE[1,i] + Flowering.DE[3,i]) / Flowering.DE[4,i] * 100),2))
}
row.names(Flowering.DE) <- c("# down regulated genes","# unaffected genes","# up regulated genes","Total expressed genes", "Total differentially expressed genes", "% differentially expressed")

Flowering.DE
### Extracting names of differentially expressed genes between genotypes

DEgene.gt.100bp.paired <- topTags(DAvsDO.100bp.paired.lrt,n = Inf)$table[topTags(DAvsDO.100bp.paired.lrt,n = Inf)$table$FDR<0.01,]
DEgene.gt.50bp.paired <- topTags(DAvsDO.50bp.paired.lrt,n = Inf)$table[topTags(DAvsDO.50bp.paired.lrt,n = Inf)$table$FDR<0.01,]
DEgene.gt.100bp.single <- topTags(DAvsDO.100bp.single.lrt,n = Inf)$table[topTags(DAvsDO.100bp.single.lrt,n = Inf)$table$FDR<0.01,]
DEgene.gt.50bp.single <- topTags(DAvsDO.50bp.single.lrt,n = Inf)$table[topTags(DAvsDO.50bp.single.lrt,n = Inf)$table$FDR<0.01,]

write.csv(DEgene.gt.50bp.paired, "DEgene.gt.50bp.paired")
write.csv(DEgene.gt.100bp.paired, "DEgene.gt.100bp.paired")
write.csv(DEgene.gt.50bp.single, "DEgene.gt.50bp.single")
write.csv(DEgene.gt.100bp.single, "DEgene.gt.100bp.single")

### Breaking DEGs from above based on up and down regulation

up.DEgene.gt.100bp.paired <- row.names(subset(DEgene.gt.100bp.paired, logFC > 0))
down.DEgene.gt.100bp.paired <- row.names(subset(DEgene.gt.100bp.paired, logFC < 0))
up.DEgene.gt.50bp.paired <- row.names(subset(DEgene.gt.50bp.paired, logFC > 0))
down.DEgene.gt.50bp.paired <- row.names(subset(DEgene.gt.50bp.paired, logFC < 0))
up.DEgene.gt.100bp.single <- row.names(subset(DEgene.gt.100bp.single, logFC > 0))
down.DEgene.gt.100bp.single <- row.names(subset(DEgene.gt.100bp.single, logFC < 0))
up.DEgene.gt.50bp.single <- row.names(subset(DEgene.gt.50bp.single, logFC > 0))
down.DEgene.gt.50bp.single <- row.names(subset(DEgene.gt.50bp.single, logFC < 0))


### Extracting names of similarally expressed genes between genotypes

SEgene.gt.100bp.paired <- row.names(topTags(DAvsDO.100bp.paired.lrt,n = Inf)$table[topTags(DAvsDO.100bp.paired.lrt,n = Inf)$table$FDR>=0.01,])
SEgene.gt.50bp.paired <- row.names(topTags(DAvsDO.50bp.paired.lrt,n = Inf)$table[topTags(DAvsDO.50bp.paired.lrt,n = Inf)$table$FDR>=0.01,])
SEgene.gt.100bp.single <- row.names(topTags(DAvsDO.100bp.single.lrt,n = Inf)$table[topTags(DAvsDO.100bp.single.lrt,n = Inf)$table$FDR>=0.01,])
SEgene.gt.50bp.single <- row.names(topTags(DAvsDO.50bp.single.lrt,n = Inf)$table[topTags(DAvsDO.50bp.single.lrt,n = Inf)$table$FDR>=0.01,])

### Comparing names of up regulated genes
#100bp paired vs 50bp paired
length(intersect(up.DEgene.gt.100bp.paired,up.DEgene.gt.50bp.paired)) #2649 Shared up regulated genes
length(setdiff(up.DEgene.gt.100bp.paired,up.DEgene.gt.50bp.paired)) #281 genes unique to 100bp paired
length(setdiff(up.DEgene.gt.50bp.paired,up.DEgene.gt.100bp.paired)) #697 genes unique to 50bp paired
#100bp paired vs 100bp single
length(intersect(up.DEgene.gt.100bp.paired,up.DEgene.gt.100bp.single)) #2721 Shared up regulated genes
length(setdiff(up.DEgene.gt.100bp.paired,up.DEgene.gt.100bp.single)) #209 genes unique to 100bp paired
length(setdiff(up.DEgene.gt.100bp.single,up.DEgene.gt.100bp.paired)) #620 genes unique to 50bp paired
#100bp single vs 50bp paired
length(intersect(up.DEgene.gt.100bp.single,up.DEgene.gt.50bp.paired)) #2907 Shared up regulated genes
length(setdiff(up.DEgene.gt.100bp.single,up.DEgene.gt.50bp.paired)) #434 genes unique to 100bp single
length(setdiff(up.DEgene.gt.50bp.paired,up.DEgene.gt.100bp.single)) #439 genes unique to 50bp paired

### Comparing names of down regulated genes
#100bp paired vs 50bp paired
length(intersect(down.DEgene.gt.100bp.paired,down.DEgene.gt.50bp.paired)) #2450 Shared down regulated genes
length(setdiff(down.DEgene.gt.100bp.paired,down.DEgene.gt.50bp.paired)) #326 genes unique to 100bp paired
length(setdiff(down.DEgene.gt.50bp.paired,down.DEgene.gt.100bp.paired)) #649 genes unique to 50bp paired
#100bp paired vs 100bp single
length(intersect(down.DEgene.gt.100bp.paired,down.DEgene.gt.100bp.single)) #2531 Shared down regulated genes
length(setdiff(down.DEgene.gt.100bp.paired,down.DEgene.gt.100bp.single)) #245 genes unique to 100bp paired
length(setdiff(down.DEgene.gt.100bp.single,down.DEgene.gt.100bp.paired)) #545 genes unique to 50bp paired
#100bp single vs 50bp paired
length(intersect(down.DEgene.gt.100bp.single,down.DEgene.gt.50bp.paired)) #2649 Shared down regulated genes
length(setdiff(down.DEgene.gt.100bp.single,down.DEgene.gt.50bp.paired)) #427 genes unique to 100bp single
length(setdiff(down.DEgene.gt.50bp.paired,down.DEgene.gt.100bp.single)) #450 genes unique to 50bp paired

### Compairing names of similarilly expressed genes
#100bp paired vs 50bp paired
length(intersect(SEgene.gt.100bp.paired,SEgene.gt.50bp.paired)) #49446 Shared similarially regulated genes
length(setdiff(SEgene.gt.100bp.paired,SEgene.gt.50bp.paired)) #1346 genes unique to 100bp paired
length(setdiff(SEgene.gt.50bp.paired,SEgene.gt.100bp.paired)) #607 genes unique to 50bp paired
#100bp paired vs 100bp single
length(intersect(SEgene.gt.100bp.paired,SEgene.gt.100bp.single)) #49627 Shared similarially regulated genes
length(setdiff(SEgene.gt.100bp.paired,SEgene.gt.100bp.single)) #1165 genes unique to 100bp paired
length(setdiff(SEgene.gt.100bp.single,SEgene.gt.100bp.paired)) #454 genes unique to 50bp paired
#100bp single vs 50bp paired
length(intersect(SEgene.gt.100bp.single,SEgene.gt.50bp.paired)) #4919 Shared similarially regulated genes
length(setdiff(SEgene.gt.100bp.single,SEgene.gt.50bp.paired)) #889 genes unique to 100bp single
length(setdiff(SEgene.gt.50bp.paired,SEgene.gt.100bp.single)) #861 genes unique to 50bp paired


### Barplots
# bring your data to long format as needed by ggplot
small.Flowering.DE <- Flowering.DE[c(1:3,5),1:4]
melt.Flowering.DE <- melt(t(small.Flowering.DE))
melt.Flowering.DE$value <- as.numeric(as.character(melt.Flowering.DE$value))
colnames(melt.Flowering.DE) <- c("Condition","Group","Count")
# plot and facet by Condition
ggplot(data=melt.Flowering.DE, aes(Condition,Count, fill=Condition)) + facet_wrap(~ Group, scales = "free_y") + geom_bar(stat = "identity")
ggsave("DE.barplot.png", width = 20, height = 8)
# log transform and replot
log2.melt.Flowering.DE <- melt.Flowering.DE
log2.melt.Flowering.DE[,3] <- log(log2.melt.Flowering.DE[3], 2)
ggplot(data=log2.melt.Flowering.DE, aes(Condition,Count, fill=Condition)) + facet_wrap(~ Group, scale scales = "free_y") + geom_bar(stat = "identity")
