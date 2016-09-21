#!/bin/bash

cd ~/KIAT/STAR/50bp/single_end/50bp-05-11-Final-8ul_1.fq.star.dir/ &&
STAR --genomeDir /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/star_genome -- readFilesIn /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/raw_data/50bp-05-11-Final-8ul_1.fq --runThreadN 6 --outReadsUnmapped Fastx &&
cd ~/KIAT/STAR/50bp/single_end/50bp-05-11-Final-8ul_2.fq.star.dir/ &&
STAR --genomeDir /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/star_genome -- readFilesIn /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/raw_data/50bp-05-11-Final-8ul_2.fq --runThreadN 6 --outReadsUnmapped Fastx &&
cd ~/KIAT/STAR/100bp/single_end/100bp-05-11-Final-8ul_1.fq.star.dir/ &&
STAR --genomeDir /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/star_genome --readFilesIn /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data/05-11-Final-8ul_1.fq --runThreadN 6 --outReadsUnmapped Fastx &&
cd ~/KIAT/STAR/100bp/single_end/100bp-05-11-Final-8ul_2.fq.star.dir/ &&
STAR --genomeDir /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/star_genome --readFilesIn /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data/05-11-Final-8ul_2.fq --runThreadN 6 --outReadsUnmapped Fastx &&
cd ~/KIAT/STAR/100bp/paired_end/100bp-05-11-Final-8ul.star.dir/ &&
STAR --genomeDir /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/star_genome --readFilesIn /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data/05-11-Final-8ul_1.fq /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data/05-11-Final-8ul_2.fq --runThreadN 6 --outReadsUnmapped Fastx &&
cd


