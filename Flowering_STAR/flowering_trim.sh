#!/bin/bash
files="Ae_Gae_2 Ae_Gae_3 All1_Gae_2 All1_Gae_3 2 6"

echo $files
for i in $files
	do
	java -jar /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 4 /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data/${i}_1.fq /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data/${i}_2.fq paired_${i}_1.fq unpaired_${i}_1.fq paired_${i}_2.fq unpaired_${i}_2.fq ILLUMINACLIP:/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data/trim_paired_end/with_Bradseq_adapter/Bradseq_adapter.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:101 CROP:50
done
