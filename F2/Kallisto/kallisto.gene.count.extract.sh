#!/bin/bash

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_fall/F2/raw_data/bioftp.org/TBD160783_20161129

Samples=`ls -d Sample_138* | sed 's/Sample_//'`

for i in $Samples
	do
	awk '{print $4}' Sample_${i}/abundance.tsv | sed 1d | sed '1 i\'${i}'' > /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/F2/Kallisto/${i}
	done

awk '{print $1}' Sample_138-1/abundance.tsv | sed 's/target_id/ID/' > /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/F2/Kallisto/138

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/F2/Kallisto

paste 138* > F2_kallisto_gene_counts.tsv
rm 138*

echo Done
