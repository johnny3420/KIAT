#!/bin/bash
reads=/Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data
files="Ae_Gae_2 Ae_Gae_3 All1_Gae_2 All1_Gae_3 2 6"

echo $reads
echo $files
cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/raw_data/Flowering
for i in $files
	do
	fastx_trimmer -Q33 -l 50 -i ${reads}/${i}_1.fq -o 50bp_${i}_1.fq
	fastx_trimmer -Q33 -l 50 -i ${reads}/${i}_2.fq -o 50bp_${i}_2.fq
done
