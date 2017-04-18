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
f2.gq <- extract.gt(f2.vcfr, element = 'GQ', IDtoRowNames = F, as.numeric = T)


# Look at GQ distribution

gq.averages <- rowMeans(f2.gq, na.rm = T)
sum(gq.averages > 30)
sum(gq.averages > 20)

# Filtering for GQ

f2.gq.20 <- f2.gq >= 20
f2.gq.30 <- f2.gq >= 30
f2.gt.20 <- f2.gt
f2.gt.30 <- f2.gt

for (i in 1:length(f2.gq.20)){
  if(f2.gq.20[i] == FALSE || is.na(f2.gq.20[i])){
    f2.gt.20[i] <- NA
  }
}

for (i in 1:length(f2.gq.30)){
  if(f2.gq.30[i] == FALSE || is.na(f2.gq.30[i])){
    f2.gt.30[i] <- NA
  }
}

#Build data frames
f2.20 <- data.frame(cbind(f2.chrom,f2.pos,f2.gt.20), stringsAsFactors = F)
colnames(f2.20)[1:2] <- c("CHROM","POS")
colnames(f2.20) <- sub("X","",colnames(f2.20))

f2.30 <- data.frame(cbind(f2.chrom,f2.pos,f2.gt.30), stringsAsFactors = F)
colnames(f2.30)[1:2] <- c("CHROM","POS")
colnames(f2.30) <- sub("X","",colnames(f2.30))

#Combine data frames and collect stats

combined.20 <- merge(parents,f2.20)
combined.20[combined.20 == "NA"] <- NA
combined.20$Like.Ae <- rowSums(combined.20 == combined.20$Ae.gt, na.rm = T) - 1
combined.20$Like.Ol <- rowSums(combined.20 == combined.20$Ol.gt, na.rm = T) - 1
combined.20$na_count <- apply(combined.20, 1, function(x) sum(is.na(x)))
combined.20$Het <- ncol(combined.20) - 9 - combined.20$Like.Ae - combined.20$Like.Ol - combined.20$na_count
combined.20$F2_SNP_Total <- rowSums(combined.20[,c("Like.Ae","Like.Ol","Het")])
combined.20$F2_SNP_Total_No_Het <- rowSums(combined.20[,c("Like.Ae","Like.Ol")])
combined.20 <- subset(combined.20, select = -na_count)

combined.30 <- merge(parents,f2.30)
combined.30[combined.30 == "NA"] <- NA
combined.30$Like.Ae <- rowSums(combined.30 == combined.30$Ae.gt, na.rm = T) - 1
combined.30$Like.Ol <- rowSums(combined.30 == combined.30$Ol.gt, na.rm = T) - 1
combined.30$na_count <- apply(combined.30, 1, function(x) sum(is.na(x)))
combined.30$Het <- ncol(combined.30) - 9 - combined.30$Like.Ae - combined.30$Like.Ol - combined.30$na_count
combined.30$F2_SNP_Total <- rowSums(combined.30[,c("Like.Ae","Like.Ol","Het")])
combined.30$F2_SNP_Total_No_Het <- rowSums(combined.30[,c("Like.Ae","Like.Ol")])
combined.30 <- subset(combined.30, select = -na_count)

## Save Data Frames

#write.csv(combined.20, "All_SNPs_GQ20.csv",row.names = F)
#write.csv(combined.30, "All_SNPs_GQ30.csv",row.names = F)

#Extract frequencies

combined.20.freq <- cbind(combined.20[,1:2],
                       transmute(combined.20, Frac.Ae.t = Like.Ae / F2_SNP_Total,
                                 Frac.Ol.t = Like.Ol / F2_SNP_Total,
                                 Frac.Het.t = Het / F2_SNP_Total,
                                 Frac.Ae.p = Like.Ae / F2_SNP_Total_No_Het,
                                 Frac.Ol.p = Like.Ol / F2_SNP_Total_No_Het,
                                 Ratio.Ae.vs.Ol.t = Frac.Ae.t / Frac.Ol.t,
                                 Ratio.Ae.vs.Ol.p = Frac.Ae.p / Frac.Ol.p))

combined.30.freq <- cbind(combined.30[,1:2],
                          transmute(combined.30, Frac.Ae.t = Like.Ae / F2_SNP_Total,
                                    Frac.Ol.t = Like.Ol / F2_SNP_Total,
                                    Frac.Het.t = Het / F2_SNP_Total,
                                    Frac.Ae.p = Like.Ae / F2_SNP_Total_No_Het,
                                    Frac.Ol.p = Like.Ol / F2_SNP_Total_No_Het,
                                    Ratio.Ae.vs.Ol.t = Frac.Ae.t / Frac.Ol.t,
                                    Ratio.Ae.vs.Ol.p = Frac.Ae.p / Frac.Ol.p))



#Plot

### Da-Ae vs Da-Ol-1 Barplot

test <- combined.20.freq[,c(1,6,7)]

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

#ggsave("F2_20_Parent_Barplot.png")

test <- combined.30.freq[,c(1,6,7)]

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

#ggsave("F2_30_Parent_Barplot.png")

### Sum Plots

Ae <- mean(combined.20.freq$Frac.Ae.t, na.rm=T) 
Ol <- mean(combined.20.freq$Frac.Ol.t, na.rm=T)
Het <- mean(combined.20.freq$Frac.Het, na.rm=T)
all <- cbind(Ae,Ol,Het)
A <- all
View(all)
barplot(all)
p.Ae <- mean(combined.20.freq$Frac.Ae.p, na.rm=T)
p.Ol <- mean(combined.20.freq$Frac.Ol.p, na.rm=T)
p <- cbind(p.Ae,p.Ol)
barplot(p)

Ae <- mean(combined.30.freq$Frac.Ae.t, na.rm=T) 
Ol <- mean(combined.30.freq$Frac.Ol.t, na.rm=T)
Het <- mean(combined.30.freq$Frac.Het, na.rm=T)
all <- cbind(Ae,Ol,Het)
B <- all
View(all)
barplot(all)
p.Ae <- mean(combined.30.freq$Frac.Ae.p, na.rm=T)
p.Ol <- mean(combined.30.freq$Frac.Ol.p, na.rm=T)
p <- cbind(p.Ae,p.Ol)
barplot(p)

### Random Stats
# Number of sites
nrow(combined.20.freq) #36149
nrow(combined.30.freq) #36149

### Filter for at least median F2 sample calls per site
median(combined.20$F2_SNP_Total) #88
median(combined.30$F2_SNP_Total) #77

combined.20 <- subset(combined.20, F2_SNP_Total >= 88)
combined.30 <- subset(combined.20, F2_SNP_Total >= 77)

# Saving combined.30
temp <- combined.30[,-c(173:177)]
write.table(temp,file = "F2_Final_SNP_Calls", col.names = T, row.names = F, sep = "\t")
test <- read.delim("F2_Final_SNP_Calls", header = T, check.names = F)
test
#Extract frequencies

combined.20.freq <- cbind(combined.20[,1:2],
                          transmute(combined.20, Frac.Ae.t = Like.Ae / F2_SNP_Total,
                                    Frac.Ol.t = Like.Ol / F2_SNP_Total,
                                    Frac.Het.t = Het / F2_SNP_Total,
                                    Frac.Ae.p = Like.Ae / F2_SNP_Total_No_Het,
                                    Frac.Ol.p = Like.Ol / F2_SNP_Total_No_Het,
                                    Ratio.Ae.vs.Ol.t = Frac.Ae.t / Frac.Ol.t,
                                    Ratio.Ae.vs.Ol.p = Frac.Ae.p / Frac.Ol.p))

combined.30.freq <- cbind(combined.30[,1:2],
                          transmute(combined.30, Frac.Ae.t = Like.Ae / F2_SNP_Total,
                                    Frac.Ol.t = Like.Ol / F2_SNP_Total,
                                    Frac.Het.t = Het / F2_SNP_Total,
                                    Frac.Ae.p = Like.Ae / F2_SNP_Total_No_Het,
                                    Frac.Ol.p = Like.Ol / F2_SNP_Total_No_Het,
                                    Ratio.Ae.vs.Ol.t = Frac.Ae.t / Frac.Ol.t,
                                    Ratio.Ae.vs.Ol.p = Frac.Ae.p / Frac.Ol.p))



#Plot

### Da-Ae vs Da-Ol-1 Barplot

test <- combined.20.freq[,c(1,6,7)]

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

#ggsave("F2_20_Parent_Barplot.png")

test <- combined.30.freq[,c(1,6,7)]

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

#ggsave("F2_30_Parent_Barplot.png")

### Sum Plots

Ae <- mean(combined.20.freq$Frac.Ae.t, na.rm=T) 
Ol <- mean(combined.20.freq$Frac.Ol.t, na.rm=T)
Het <- mean(combined.20.freq$Frac.Het, na.rm=T)
all <- cbind(Ae,Ol,Het)
C <- all
View(all)
barplot(all)
p.Ae <- mean(combined.20.freq$Frac.Ae.p, na.rm=T)
p.Ol <- mean(combined.20.freq$Frac.Ol.p, na.rm=T)
p <- cbind(p.Ae,p.Ol)
barplot(p)

Ae <- mean(combined.30.freq$Frac.Ae.t, na.rm=T) 
Ol <- mean(combined.30.freq$Frac.Ol.t, na.rm=T)
Het <- mean(combined.30.freq$Frac.Het, na.rm=T)
all <- cbind(Ae,Ol,Het)
D <- all
View(all)
barplot(all)
p.Ae <- mean(combined.30.freq$Frac.Ae.p, na.rm=T)
p.Ol <- mean(combined.30.freq$Frac.Ol.p, na.rm=T)
p <- cbind(p.Ae,p.Ol)
barplot(p)

### Combined plots
Ratios <- rbind(A,B,C,D)
rownames(Ratios) <- c("GQ_20","GQ_30","GQ_20_Median_Filtered","GQ_30_Median_Filtered")
m.Ratios <- melt(Ratios)
colnames(m.Ratios) <- c("Filter","Genotype","Proportion")
pl <- ggplot(m.Ratios, aes(Genotype,Proportion,fill = Genotype))
pl + geom_col() + facet_grid(~Filter)
#ggsave("Filter_Comparision.png")
### Extract SNP Postions

positions <- combined.30[,1:2] #18226
#write.csv(positions, "Final_F2_SNP_Sites.csv", row.names = F)
#write.table(positions, "Final_F2_SNP_Sites.tab", row.names = F, sep = "\t")

### Napus Rna Seq dataset avail, check to see what they found

### Look at SNPs Ruijuan has but i dont in IGV