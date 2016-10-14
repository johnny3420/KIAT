#!/bin/bash

files="Ae_Gae_2 Ae_Gae_3 All1_Gae_2 All1_Gae_3 2 6"

for i in $files
	do
	awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data/${i}_paired_1.fq > length.check/${i}_paired_1.length
	awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data/${i}_paired_2.fq > length.check/${i}_paired_2.length
	done
echo DoNe
