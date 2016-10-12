#!/bin/bash
files="Ae_Gae_2 Ae_Gae_3 All1_Gae_2 All1_Gae_3 2 6"
echo $files
cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/Flowering_STAR/Alignments
for i in $files
	do
	if [[ -d short.${i}.paired.star.dir ]]
		then
			echo OK
	else
		mkdir short.${i}.paired.star.dir
	fi

	cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/Flowering_STAR/Alignments/short.${i}.paired.star.dir
	STAR --genomeDir /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/star_genome --readFilesIn ../../Trimmed_Reads/${i}_1.fq ../../Trimmed_Reads/${i}_2.fq --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus.annotation_v5_modified_modified.gff3 --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic –alignIntronMax 15000 --outFilterIntronMotifs RemoveNoncanonical --runThreadN 5 --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfeatureExon CDS --outReadsUnmapped Fastx
	cd ..
	
	if [[ -d short.${i}.single.star.dir ]]
		then
			echo OK
	else
		mkdir short.${i}.single.star.dir
	fi
	cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/Flowering_STAR/Alignments/short.${i}.single.star.dir
	STAR --genomeDir /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/star_genome --readFilesIn ../../Trimmed_Reads/${i}_1.fq --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus.annotation_v5_modified_modified.gff3 --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic –alignIntronMax 15000 --outFilterIntronMotifs RemoveNoncanonical --runThreadN 5 --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfeatureExon CDS --outReadsUnmapped Fastx
	cd ..

done

echo DoNe
		
