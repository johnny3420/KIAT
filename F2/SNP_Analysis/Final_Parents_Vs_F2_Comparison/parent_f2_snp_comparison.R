### F2 and Parent SNP comparison

library(vcfR)
library(dplyr)
library(ggplot2)
library(reshape2)

#read in and format parent data

parents <- read.csv("vcf.Ae.Ol.intersect.final.csv", stringsAsFactors = F)
parents <- parents[,-1]

#remove 'random' chromosomes

random <- unique(grep("_random", parents$CHROM, value = T))
parents <- parents[!(parents$CHROM %in% random),]

#read in and format f2 data

f2.vcfr <- read.vcfR("F2_filtered.vcf.gz")
f2.chrom <- getCHROM(f2.vcfr)
f2.pos <- getPOS(f2.vcfr)
f2.gt <- extract.gt(f2.vcfr, element = 'GT', IDtoRowNames = F)
f2.gl <- extract.gt(f2.vcfr, element = 'GL', IDtoRowNames = F)


## Filtering 
# filter for GL

for(i in 1:length(f2.gl)){
  if(!is.na(f2.gl[i])){
    if(sum(sapply(as.numeric(unlist(strsplit(f2.gl[i],","))), function(x) (x > -2 && x != 0))) >= 1){
      f2.gl[i] <- NA
    }
  }
}

for(i in 1:length(f2.gl)){
  if(is.na(f2.gl[i])){
    f2.gt[i] <- NA
  }
}



f2 <- data.frame(cbind(f2.chrom,f2.pos,f2.gt), stringsAsFactors = F)
colnames(f2)[1:2] <- c("CHROM","POS")
colnames(f2) <- sub("X","",colnames(f2))

#Combine data frames and collect stats

combined <- merge(parents,f2)
combined[combined == "NA"] <- NA
combined$Like.Ae <- rowSums(combined == combined$Ae.gt, na.rm = T) - 1
combined$Like.Ol <- rowSums(combined == combined$Ol.gt, na.rm = T) - 1
combined$na_count <- apply(combined, 1, function(x) sum(is.na(x)))
combined$Het <- ncol(combined) - 9 - combined$Like.Ae - combined$Like.Ol - combined$na_count
combined$F2_SNP_Total <- rowSums(combined[,c("Like.Ae","Like.Ol","Het")])
combined$F2_SNP_Total_No_Het <- rowSums(combined[,c("Like.Ae","Like.Ol")])
combined <- subset(combined, select = -na_count)

### Filter for at least 50 F2 sample calls per site

filtered.combined <- subset(combined, F2_SNP_Total >= 43)

#Extract frequencies

combined.freq <- cbind(combined[,1:2],
                       transmute(combined, Frac.Ae.t = Like.Ae / F2_SNP_Total,
                                 Frac.Ol.t = Like.Ol / F2_SNP_Total,
                                 Frac.Het.t = Het / F2_SNP_Total,
                                 Frac.Ae.p = Like.Ae / F2_SNP_Total_No_Het,
                                 Frac.Ol.p = Like.Ol / F2_SNP_Total_No_Het,
                                 Ratio.Ae.vs.Ol.t = Frac.Ae.t / Frac.Ol.t,
                                 Ratio.Ae.vs.Ol.p = Frac.Ae.p / Frac.Ol.p))



#Plot

### Da-Ae vs Da-Ol-1 Barplot

test <- combined.freq[,c(1,6,7)]

test2 <- aggregate(.~CHROM, data=test, mean)
test3 <- c("Overall",mean(test2$Frac.Ae.p),mean(test2$Frac.Ol.p))
test2 <- data.frame(rbind(test2,test3),stringsAsFactors=F)
colnames(test2)[2:3] <- c("Da-Ae", "Da-Ol-1")
test2$`Da-Ae` <- as.numeric(as.character(test2$`Da-Ae`))
test2$`Da-Ol-1` <- as.numeric(as.character(test2$`Da-Ol-1`))


m.test2 <- melt(test2)
colnames(m.test2)[2] <- "Genotype"

ggplot(data = m.test2, aes(x=CHROM,y=value)) +
  geom_bar(aes(fill = Genotype),stat = "identity",position = "dodge") +
  geom_hline(yintercept= .5) +
  scale_y_continuous(limits = c(0,.6), breaks = c(.1,.2,.3,.4,.5,.6)) +
  labs(list(x = "Chromosome", y = "Proportion"))

ggsave("F2_Parent_Barplot.png")

### Sum Plots

Ae <- mean(combined.freq$Frac.Ae.t, na.rm=T) 
Ol <- mean(combined.freq$Frac.Ol.t, na.rm=T)
Het <- mean(combined.freq$Frac.Het, na.rm=T)
all <- cbind(Ae,Ol,Het)
View(all)
barplot(all)
p.Ae <- mean(combined.freq$Frac.Ae.p, na.rm=T)
p.Ol <- mean(combined.freq$Frac.Ol.p, na.rm=T)
p <- cbind(p.Ae,p.Ol)
barplot(p)

### Napus Rna Seq dataset avail, check to see what they found

### Look at SNPs Ruijuan has but i dont in IGV