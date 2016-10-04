#!/bin/bash

bwa mem -t 4 Brassica_napus_v4.1.chromosomes.fa ../raw_data/Flowering/100bp_2_1.fq ../raw_data/Flowering/100bp_2_2.fq > 2_100bp.paired/aln.sam &&
bwa mem -t 4 Brassica_napus_v4.1.chromosomes.fa ../raw_data/Flowering/100bp_2_1.fq > 2_100bp.single/aln.sam &&
bwa mem Brassica_napus_v4.1.chromosomes.fa ../raw_data/Flowering/50bp_2_1.fq > mem.2_50bp.paired/aln.sam &&
bwa mem Brassica_napus_v4.1.chromosomes.fa ../raw_data/Flowering/50bp_2_2.fq > mem.2_50bp.single/aln.sam &&
echo DoNe
