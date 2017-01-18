#!/bin/bash

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_fall/F2/raw_data/bioftp.org/TBD160783_20161129
Folder1=`ls -d Sample* | sed 's/S.*-//'`
echo $Folder1

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/F2/STAR_Runs

for i in $Folder1
	do
        if [[ -d Sample_138-${i}.paired.star.dir ]]
               then
                       echo OK
	else
                mkdir Sample_138-${i}.paired.star.dir
        fi

	cd Sample_138-${i}.paired.star.dir
	STAR --genomeDir /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/star_genome --readFilesIn /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_fall/F2/raw_data/bioftp.org/TBD160783_20161129/Sample_138-${i}/138-${i}_paired_1.fq.gz /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_fall/F2/raw_data/bioftp.org/TBD160783_20161129/Sample_138-${i}/138-${i}_paired_2.fq.gz --sjdbGTFfile /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus/Brassica_napus.annotation_v5_modified_modified.gff3 --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic --alignIntronMax 15000 --outFilterIntronMotifs RemoveNoncanonical --runThreadN 5 --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfeatureExon CDS --outReadsUnmapped Fastx --readFilesCommand zcat
	cd ..
	echo Sample_138-${i} is done
	done

echo DOONE
