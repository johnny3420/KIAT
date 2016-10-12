#!/bin/bash
cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/Flowering_Star/Alignments
files="Ae_Gae_2 Ae_Gae_3 All1_Gae_2 All1_Gae_3 2 6"
echo $files
for i in $files
	do
	cd long.${i}.paired.star.dir
	awk '{print $1 "\t" $2}' ReadsPerGene.out.tab| sed 1,4d | sed '1 i\ID '100bp_paired_${i}'' > ../../Gene_Counts/100bp.paired.${i}.counts
	cd ../long.${i}.single.star.dir
        awk '{print $1 "\t" $2}' ReadsPerGene.out.tab| sed 1,4d | sed '1 i\ID '100bp_single_${i}'' > ../../Gene_Counts/100bp.single.${i}.counts
	cd ../short.${i}.paired.star.dir
	awk '{print $1 "\t" $2}' ReadsPerGene.out.tab| sed 1,4d | sed '1 i\ID '50bp_paired_${i}'' > ../../Gene_Counts/50bp.paired.${i}.counts
	cd ../short.${i}.single.star.dir
	awk '{print $1 "\t" $2}' ReadsPerGene.out.tab| sed 1,4d | sed '1 i\ID '50bp_single_${i}'' > ../../Gene_Counts/50bp.single.${i}.counts
	cd ..	
	done
cd ../Gene_Counts
paste 50bp* | awk '{ for (i=1;i<=NF;i+=2) {$i=""} print $0}' | sed 's/ /\t/g' > 50bp.reads
paste 50bp* | awk '{print $1}' > 50bp.ID
paste 50bp.ID 50bp.reads > flowering.star.50bp.read.count.tsv
paste 100bp* | awk '{ for (i=1;i<=NF;i+=2) {$i=""} print $0}' | sed 's/ /\t/g' > 100bp.reads
paste 100bp* | awk '{print $1}' > 100bp.ID
paste 100bp.reads > flowering.star.100bp.read.count.tsv
paste flowering* > flowering.star.read.count.tsv

