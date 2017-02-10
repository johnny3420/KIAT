#!/bin/bash
files="Ae_Gae_2 Ae_Gae_3 All1_Gae_2 All1_Gae_3 2 6"
echo $files
cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/Flowering_STAR/Alignments/100bp_Alignments
for i in $files
	do
        if [[ -d long.${i}.paired.star.dir ]]
                then
                        echo OK
        else
                mkdir long.${i}.paired.star.dir
        fi
	cp -r /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data/${i}_paired.star.trim.dir/* /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/Flowering_STAR/Alignments/long.${i}.paired.star.dir
	if [[ -d long.${i}.single.star.dir ]]
		then
			echo OK
	else
		mkdir long.${i}.single.star.dir
	fi
	cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/Flowering_STAR/Alignments/long.${i}.single.star.dir
	STAR --genomeDir /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/star_genome --readFilesIn /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data/${i}_paired_1.fq --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus.annotation_v5_modified_modified.gff3 --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic --alignIntronMax 15000 --outFilterIntronMotifs RemoveNoncanonical --runThreadN 5 --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfeatureExon CDS --outReadsUnmapped Fastx
	cd ..

done

echo DoNe
		
