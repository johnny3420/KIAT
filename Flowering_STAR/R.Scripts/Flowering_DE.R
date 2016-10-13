### Comparing different conditions

#Load files and library
library(edgeR)
library(reshape2)
library(ggplot2)
library(plyr)
library(gplots)
load("R_Analysis/flowering.read.count.small.Rdata")
load("R_Analysis/flowering.read.count.sample.Rdata")

# Create DGElist and Normalize
dge.data.flowering <- DGEList(counts=flowering.read.count.small, group=flowering.read.count.sample$group)
dge.data.flowering <- calcNormFactors(dge.data.flowering, method = "TMM") 

# Setting up design
design.flowering <- model.matrix(~0+group, data=dge.data.flowering$samples)
colnames(design.flowering) <- levels(dge.data.flowering$samples$group)
design.flowering

# estimating dispersions and create fit
dge.data.flowering <- estimateGLMCommonDisp(dge.data.flowering, design.flowering,verbose = TRUE) # Disp = 0.15732 , BCV = 0.3966
dge.data.flowering <- estimateGLMTrendedDisp(dge.data.flowering,design.flowering)
dge.data.flowering <- estimateGLMTagwiseDisp(dge.data.flowering,design.flowering)
plotBCV(dge.data.flowering)
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
row.names(Flowering.DE) <- c("# down regulated genes","# No differentiation","# up regulated genes","Total expressed genes", "Total differentially expressed genes", "% differentially expressed")

Flowering.DE
write.csv(Flowering.DE, "R_Analysis/Flowering.DE.csv")
### Extracting differentially expressed genes between genotypes

DEgene.gt.100bp.paired <- topTags(DAvsDO.100bp.paired.lrt,n = Inf)$table[topTags(DAvsDO.100bp.paired.lrt,n = Inf)$table$FDR<0.01,]
DEgene.gt.50bp.paired <- topTags(DAvsDO.50bp.paired.lrt,n = Inf)$table[topTags(DAvsDO.50bp.paired.lrt,n = Inf)$table$FDR<0.01,]
DEgene.gt.100bp.single <- topTags(DAvsDO.100bp.single.lrt,n = Inf)$table[topTags(DAvsDO.100bp.single.lrt,n = Inf)$table$FDR<0.01,]
DEgene.gt.50bp.single <- topTags(DAvsDO.50bp.single.lrt,n = Inf)$table[topTags(DAvsDO.50bp.single.lrt,n = Inf)$table$FDR<0.01,]

write.csv(DEgene.gt.50bp.paired, "R_Analysis/DEgene.gt.50bp.paired.csv")
write.csv(DEgene.gt.100bp.paired, "R_Analysis/DEgene.gt.100bp.paired.csv")
write.csv(DEgene.gt.50bp.single, "R_Analysis/DEgene.gt.50bp.single.csv")
write.csv(DEgene.gt.100bp.single, "R_Analysis/DEgene.gt.100bp.single.csv")

###Combining into 1 tab delim file of only gene names
p100 <- as.data.frame(row.names(DEgene.gt.100bp.paired))
p50 <- as.data.frame(row.names(DEgene.gt.50bp.paired))
s100 <- as.data.frame(row.names(DEgene.gt.100bp.single))
s50 <- as.data.frame(row.names(DEgene.gt.50bp.single))
cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
  rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}
combined <- data.frame(cbind.fill(p100,p50,s100,s50))
colnames(combined) <- c("100bp_paired","50bp_paired","100bp_single","50bp_single")
write.table(combined, "R_Analysis/Combined_DE_gene_names.tab")

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
length(intersect(up.DEgene.gt.100bp.paired,up.DEgene.gt.50bp.paired)) #2554 Shared up regulated genes
length(setdiff(up.DEgene.gt.100bp.paired,up.DEgene.gt.50bp.paired)) #250 genes unique to 100bp paired
length(setdiff(up.DEgene.gt.50bp.paired,up.DEgene.gt.100bp.paired)) #504 genes unique to 50bp paired
#100bp paired vs 100bp single
length(intersect(up.DEgene.gt.100bp.paired,up.DEgene.gt.100bp.single)) #2589 Shared up regulated genes
length(setdiff(up.DEgene.gt.100bp.paired,up.DEgene.gt.100bp.single)) #215 genes unique to 100bp paired
length(setdiff(up.DEgene.gt.100bp.single,up.DEgene.gt.100bp.paired)) #452 genes unique to 50bp paired
#100bp single vs 50bp paired
length(intersect(up.DEgene.gt.100bp.single,up.DEgene.gt.50bp.paired)) #2690 Shared up regulated genes
length(setdiff(up.DEgene.gt.100bp.single,up.DEgene.gt.50bp.paired)) #351 genes unique to 100bp single
length(setdiff(up.DEgene.gt.50bp.paired,up.DEgene.gt.100bp.single)) #368 genes unique to 50bp paired
#100bp paired vs 50bp single
length(intersect(up.DEgene.gt.100bp.paired,up.DEgene.gt.50bp.single)) #2235 Shared up regulated genes
length(setdiff(up.DEgene.gt.100bp.paired,up.DEgene.gt.50bp.single)) #569 genes unique to 100bp paired
length(setdiff(up.DEgene.gt.50bp.single,up.DEgene.gt.100bp.paired)) #533 genes unique to 50bp single
#100bp single vs 50bp single
length(intersect(up.DEgene.gt.100bp.single,up.DEgene.gt.50bp.single)) #2427 Shared up regulated genes
length(setdiff(up.DEgene.gt.100bp.single,up.DEgene.gt.50bp.single)) #614 genes unique to 100bp paired
length(setdiff(up.DEgene.gt.50bp.single,up.DEgene.gt.100bp.single)) #341 genes unique to 50bp single
#Put into a table
Only_1st <- c(250,215,351,569,614)
Only_2nd <- c(504,452,368,533,341)
Shared <- c(2554,2589,2690,2235,2427)
up.reg <- t(data.frame(Only_1st,Only_2nd,Shared))
colnames(up.reg) <-c("100p_vs_50p","100p_vs_100s","100s_vs_50p","100p_vs_50s","100s_vs_50s")

### Comparing names of down regulated genes
#100bp paired vs 50bp paired
length(intersect(down.DEgene.gt.100bp.paired,down.DEgene.gt.50bp.paired)) #2433 Shared down regulated genes
length(setdiff(down.DEgene.gt.100bp.paired,down.DEgene.gt.50bp.paired)) #253 genes unique to 100bp paired
length(setdiff(down.DEgene.gt.50bp.paired,down.DEgene.gt.100bp.paired)) #497 genes unique to 50bp paired
#100bp paired vs 100bp single
length(intersect(down.DEgene.gt.100bp.paired,down.DEgene.gt.100bp.single)) #2458 Shared down regulated genes
length(setdiff(down.DEgene.gt.100bp.paired,down.DEgene.gt.100bp.single)) #228 genes unique to 100bp paired
length(setdiff(down.DEgene.gt.100bp.single,down.DEgene.gt.100bp.paired)) #441 genes unique to 50bp paired
#100bp single vs 50bp paired
length(intersect(down.DEgene.gt.100bp.single,down.DEgene.gt.50bp.paired)) #2546 Shared down regulated genes
length(setdiff(down.DEgene.gt.100bp.single,down.DEgene.gt.50bp.paired)) #353 genes unique to 100bp single
length(setdiff(down.DEgene.gt.50bp.paired,down.DEgene.gt.100bp.single)) #384 genes unique to 50bp paired
#100bp paired vs 50bp single
length(intersect(down.DEgene.gt.100bp.paired,down.DEgene.gt.50bp.single)) #2194 Shared down regulated genes
length(setdiff(down.DEgene.gt.100bp.paired,down.DEgene.gt.50bp.single)) #492 genes unique to 100bp paired
length(setdiff(down.DEgene.gt.50bp.single,down.DEgene.gt.100bp.paired)) #608 genes unique to 50bp single
#100bp single vs 50bp single
length(intersect(down.DEgene.gt.100bp.single,down.DEgene.gt.50bp.single)) #2392 Shared down regulated genes
length(setdiff(down.DEgene.gt.100bp.single,down.DEgene.gt.50bp.single)) #507 genes unique to 100bp paired
length(setdiff(down.DEgene.gt.50bp.single,down.DEgene.gt.100bp.single)) #410 genes unique to 50bp single
# Put into table
Only_1st <- c(253,228,353,492,507)
Only_2nd <- c(497,441,384,608,410)
Shared <- c(2433,2458,2546,2194,2392)
down.reg <- t(data.frame(Only_1st,Only_2nd,Shared))
colnames(down.reg) <-c("100p_vs_50p","100p_vs_100s","100s_vs_50p","100p_vs_50s","100s_vs_50s")

## Combining all DE genes
de.reg <- up.reg + down.reg

## Plotting of DE genes between different comparisons
melt.de.reg <- melt(de.reg, varnames = c("Ownership","Comparison"), value.name = "Gene_Counts")
ggplot(data=melt.de.reg, aes(Ownership,Gene_Counts, fill=Ownership)) + ggtitle("Differently Expressed Gene Counts Between Conditions") +
  facet_wrap(~ Comparison, scales = "free_y") + geom_bar(stat = "identity")
ggsave("R_Analysis/Combined.DE.barplot.png")
### Compairing names of similarilly expressed genes
#100bp paired vs 50bp paired
length(intersect(SEgene.gt.100bp.paired,SEgene.gt.50bp.paired)) #47408 Shared similarially regulated genes
length(setdiff(SEgene.gt.100bp.paired,SEgene.gt.50bp.paired)) #1001 genes unique to 100bp paired
length(setdiff(SEgene.gt.50bp.paired,SEgene.gt.100bp.paired)) #503 genes unique to 50bp paired
#100bp paired vs 100bp single
length(intersect(SEgene.gt.100bp.paired,SEgene.gt.100bp.single)) #47516 Shared similarially regulated genes
length(setdiff(SEgene.gt.100bp.paired,SEgene.gt.100bp.single)) #893 genes unique to 100bp paired
length(setdiff(SEgene.gt.100bp.single,SEgene.gt.100bp.paired)) #443 genes unique to 50bp paired
#100bp single vs 50bp paired
length(intersect(SEgene.gt.100bp.single,SEgene.gt.50bp.paired)) #47207 Shared similarially regulated genes
length(setdiff(SEgene.gt.100bp.single,SEgene.gt.50bp.paired)) #752 genes unique to 100bp single
length(setdiff(SEgene.gt.50bp.paired,SEgene.gt.100bp.single)) #704 genes unique to 50bp paired
#100bp paired vs 50bp single
length(intersect(SEgene.gt.100bp.paired,SEgene.gt.50bp.single)) #47268 Shared similarially regulated genes
length(setdiff(SEgene.gt.100bp.paired,SEgene.gt.50bp.single)) #1141 genes unique to 100bp paired
length(setdiff(SEgene.gt.50bp.single,SEgene.gt.100bp.paired)) #1061 genes unique to 50bp single
#100bp single vs 50bp single
length(intersect(SEgene.gt.100bp.single,SEgene.gt.50bp.single)) #47208 Shared down regulated genes
length(setdiff(SEgene.gt.100bp.single,SEgene.gt.50bp.single)) #751 genes unique to 100bp paired
length(setdiff(SEgene.gt.50bp.single,SEgene.gt.100bp.single)) #1121 genes unique to 50bp single
# Put into table
Only_1st <- c(1001,893,752,1141,751)
Only_2nd <- c(503,443,704,1061,1121)
Shared <- c(47408,47516,47207,47268,47208)
sim.reg <- t(data.frame(Only_1st,Only_2nd,Shared))
colnames(sim.reg) <-c("100p_vs_50p","100p_vs_100s","100s_vs_50p","100p_vs_50s","100s_vs_50s")

### Barplots
# bring your data to long format as needed by ggplot
small.Flowering.DE <- Flowering.DE[c(1:3,5),1:4]
melt.Flowering.DE <- melt(t(small.Flowering.DE))
melt.Flowering.DE$value <- as.numeric(as.character(melt.Flowering.DE$value))
colnames(melt.Flowering.DE) <- c("Comparison","Group","Count")
melt.Flowering.DE$Comparison <- gsub("\\.", " ", melt.Flowering.DE$Comparison)
melt.Flowering.DE$Comparison <- sub("DAvsDO", "DA vs DO", melt.Flowering.DE$Comparison)
# plot and facet by Comparison
ggplot(data=melt.Flowering.DE, aes(Comparison,Count, fill=Comparison)) + facet_wrap(~ Group, scales = "free_y") + 
  geom_bar(stat = "identity") + ggtitle("Differentially Expressed Gene Counts")
ggsave("R_Analysis/DE.barplot.png", width = 20, height = 8)

### Venn Diagrams

# Up Regulated Genes
up.venn <- venn(list("100bp Paired"=up.DEgene.gt.100bp.paired,"100bp Single"=up.DEgene.gt.100bp.single,"50bp Paired"=up.DEgene.gt.50bp.paired,"50bp Single"=up.DEgene.gt.50bp.single))
# Down Regulated Genes
venn(list("100bp Paired"=down.DEgene.gt.100bp.paired,"100bp Single"=down.DEgene.gt.100bp.single,"50bp Paired"=down.DEgene.gt.50bp.paired,"50bp Single"=down.DEgene.gt.50bp.single))
# No Differentiation
venn(list("100bp Paired"=SEgene.gt.100bp.paired,"100bp Single"=SEgene.gt.100bp.single,"50bp Paired"=SEgene.gt.50bp.paired,"50bp Single"=SEgene.gt.50bp.single))

