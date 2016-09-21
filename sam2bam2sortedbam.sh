#!/bin/bash

cd ~/KIAT/STAR/50bp/single_end/50bp-05-11-Final-8ul_1.fq.star.dir/ &&
samtools view -b -S Aligned.out.sam > Aligned.out.bam &&
samtools sort Aligned.out.bam -o Aligned.out.sorted.bam &&
pwd &&
echo done &&
cd ~/KIAT/STAR/50bp/single_end/50bp-05-11-Final-8ul_2.fq.star.dir/ &&
samtools view -b -S Aligned.out.sam > Aligned.out.bam &&
samtools sort Aligned.out.bam -o Aligned.out.sorted.bam &&
pwd &&
echo done &&
cd ~/KIAT/STAR/50bp/paired_end/50bp-05-11-Final-8ul.star.dir/ &&
samtools view -b -S Aligned.out.sam > Aligned.out.bam &&
samtools sort Aligned.out.bam -o Aligned.out.sorted.bam &&
pwd &&
echo done &&
cd ~/KIAT/STAR/100bp/single_end/100bp-05-11-Final-8ul_1.fq.star.dir/ &&
samtools view -b -S Aligned.out.sam > Aligned.out.bam &&
samtools sort Aligned.out.bam -o Aligned.out.sorted.bam &&
pwd &&
echo done &&
cd ~/KIAT/STAR/100bp/single_end/100bp-05-11-Final-8ul_2.fq.star.dir/ &&
samtools view -b -S Aligned.out.sam > Aligned.out.bam &&
samtools sort Aligned.out.bam -o Aligned.out.sorted.bam &&
pwd &&
echo done &&
cd ~/KIAT/STAR/100bp/paired_end/100bp-05-11-Final-8ul.star.dir/ &&
samtools view -b -S Aligned.out.sam > Aligned.out.bam &&
samtools sort Aligned.out.bam -o Aligned.out.sorted.bam &&
pwd &&
echo done &&
cd


