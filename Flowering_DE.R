### Comparing different conditions
#Load files
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

## Finding differentially expressed genes with genotype and constant pairing, based on read length
# Da-Ae paired
DA.P.length.lrt <-glmLRT(fit, contrast = c(1,0,0,0,-1,0,0,0))
DA.P.length <- summary(decideTestsDGE(DA.P.length.lrt,p.value=0.01))
# Da-Ol-1 paired
DO.P.length.lrt <-glmLRT(fit, contrast = c(0,1,0,0,0,-1,0,0))
DO.P.length <- summary(decideTestsDGE(DO.P.length.lrt,p.value=0.01))
# Da-Ae single
DA.S.length.lrt <-glmLRT(fit, contrast = c(0,0,1,0,0,0,-1,0))
DA.S.length <- summary(decideTestsDGE(DA.S.length.lrt,p.value=0.01))
# Da-Ol-1 single
DO.S.length.lrt <-glmLRT(fit, contrast = c(0,0,0,1,0,0,0,-1))
DO.S.length <- summary(decideTestsDGE(DO.S.length.lrt,p.value=0.01))

length <- cbind(DA.P.length,DO.P.length,DA.S.length,DO.S.length)
colnames(length) <- c("DA.P.length","DO.P.length","DA.S.length","DO.S.length")

## Finding differentially expressed genes with genotype and constant length, based on pairing

# Da-Ae 100bp
DA.100bp.pairing.lrt <-glmLRT(fit, contrast = c(1,0,-1,0,0,0,0,0))
DA.100bp.pairing <- summary(decideTestsDGE(DA.P.length.lrt,p.value=0.01))
# Da-Ol-1 100bp
DO.100bp.pairing.lrt <-glmLRT(fit, contrast = c(0,1,0,-1,0,0,0,0))
DO.100bp.pairing <- summary(decideTestsDGE(DO.P.length.lrt,p.value=0.01))
# Da-Ae 50bp
DA.50bp.pairing.lrt <-glmLRT(fit, contrast = c(0,0,0,0,1,0,-1,0))
DA.50bp.pairing <- summary(decideTestsDGE(DA.S.length.lrt,p.value=0.01))
# Da-Ol-1 50bp
DO.50bp.pairing.lrt <-glmLRT(fit, contrast = c(0,0,0,0,0,1,0,-1))
DO.50bp.pairing <- summary(decideTestsDGE(DO.S.length.lrt,p.value=0.01))

pairing <- cbind(DA.100bp.pairing,DO.100bp.pairing,DA.50bp.pairing,DO.50bp.pairing)
colnames(pairing) <- c("DA.100bp.pairing","DO.100bp.pairing","DA.50bp.pairing","DO.50bp.pairing")

## Finding differentially expressed genes between different genotypes, constant pairing and length

# 100bp paired
DAvsDO.100bp.paired.lrt <-glmLRT(fit, contrast = c(
  DAvsDO.100bp.paired
  # 100bp single
  DAvsDO.100bp.single.lrt <-glmLRT(fit, contrast = c(
    DAvsDO.100bp.single
    # 50bp paired
    DAvsDO.50bp.paired.lrt <-glmLRT(fit, contrast = c(
      DAvsDO.50bp.paired
      # 50bp single
      DAvsDO.50bp.single.lrt <-glmLRT(fit, contrast = c(
        DAvsDO.50bp.single