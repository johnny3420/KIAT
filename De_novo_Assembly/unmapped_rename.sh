#!/bin/bash 

# 1) reformat unmapped & raw read 2 into 4 line format
# 2) compare two files, if the 1st & 3rd columns are the same between the two reformatted files, replace the 2nd column of unmapped file with the 2nd column of raw read file

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data

sample=`ls *_paired_1.fq | sed 's/.fq//g' | sed 's/_paired_1$//g'`

echo $sample

for i in $sample
	do
	echo ${i}
	sed 'N;s/\n/ /g' /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/De_novo_Assembly/Alignments/Mate_1/${i}.Unmapped.out.mate1 | sed 'N;s/\n/ /g' > /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/De_novo_Assembly/Alignments/Mate_1/${i}.Unmapped.out.mate1.combined
	sed 'N;s/\n/ /g' /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data/${i}_1.fq | sed 'N;s/\n/ /g' > /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/De_novo_Assembly/Alignments/Mate_1/${i}_1.fq.combined
	awk 'NR==FNR{a[$1]=$0}NR>FNR{if ($1 in a) print $0}' /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/De_novo_Assembly/Alignments/Mate_1/${i}.Unmapped.out.mate1.combined /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/De_novo_Assembly/Alignments/Mate_1/${i}_1.fq.combined | awk '{print $1, $2 "\n" $3 "\n" $4 "\n" $5}' > /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/De_novo_Assembly/Alignments/Renamed_Mate_1/${i}.Unmapped.out.mate1.replaceID
	sed 'N;s/\n/ /g' /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/De_novo_Assembly/Alignments/Mate_2/${i}.Unmapped.out.mate2 | sed 'N;s/\n/ /g' > /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/De_novo_Assembly/Alignments/Mate_2/${i}.Unmapped.out.mate2.combined
	sed 'N;s/\n/ /g' /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data/${i}_2.fq | sed 'N;s/\n/ /g' > /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/De_novo_Assembly/Alignments/Mate_2/${i}_2.fq.combined
	awk 'NR==FNR{a[$1]=$0}NR>FNR{if ($1 in a) print $0}' /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/De_novo_Assembly/Alignments/Mate_2/${i}.Unmapped.out.mate2.combined /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/De_novo_Assembly/Alignments/Mate_2/${i}_2.fq.combined | awk '{print $1, $2 "\n" $3 "\n" $4 "\n" $5}' > /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/De_novo_Assembly/Alignments/Renamed_Mate_2/${i}.Unmapped.out.mate2.replaceID
done

echo DOOONE
