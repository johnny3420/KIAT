#!/bin/bash

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data

sample=`ls *_paired_1.fq | sed 's/.fq//g' | sed 's/_1$//g'`
echo $sample

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/De_novo_Assembly/Alignments

for i in $sample
do 

	if [[ -d $i.star.trim.dir ]]
		then
			echo OK
	else 
		mkdir $i.star.trim.dir
	fi

	cd $i.star.trim.dir/

	STAR --genomeDir /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/star_genome --readFilesIn /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data/${i}_1.fq /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data/${i}_2.fq --outSAMtype BAM Unsorted --sjdbGTFfile /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus.annotation_v5_modified_modified.gff3 --twopassMode Basic --alignIntronMax 15000 --outFilterIntronMotifs RemoveNoncanonical --runThreadN 5 --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfeatureExon CDS --outReadsUnmapped Fastx
	cd ..
done
echo DOOOOOOOOOOOONE
