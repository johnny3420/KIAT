### Comparing different conditions

#Load files and library
library(edgeR)
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
  Flowering.DE[5,i] <- sprintf("%.2f", round(((Flowering.DE[1,i] + Flowering.DE[3,i]) / Flowering.DE[4,i] * 100),2))
}
row.names(Flowering.DE) <- c("# down regulated genes","# unaffected genes","# up regulated genes","Total expressed genes", "% differentially expressed")

