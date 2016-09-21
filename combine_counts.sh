#!/bin/bash

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/Kallisto
folder=`ls`
echo $folder
mkdir Read.counts
for i in $folder
	do
	cd $i
	awk '{print $1 "\t" $4}' abundance.tsv | sed 1d | sed '1 i\ID	'$i'' >> ../Read.counts/$i.abundance.tsv
	awk '{print $1 "\t" $5}' abundance.tsv | sed 1d | sed '1 i\ID	'$i'' >> ../Read.counts/$i.abundance.normalized.tsv
	cd ..
	done
cd Read.counts
paste *abundance.tsv | awk '{ for (i=1;i<=NF;i+=2) {$i=""} print $0}' | sed 's/ /\t/g' > reads
paste *abundance.normalized.tsv | awk '{ for (i=1;i<=NF;i+=2) {$i=""} print $0}' | sed 's/ /\t/g' > reads.normalized       
paste *abundance.tsv | awk '{print $1}' > ID
paste *abundance.normalized.tsv | awk '{print $1}' > ID.normalized
paste ID reads > kallisto.read.count.tsv
paste ID.normalized reads.normalized > kallisto.read.count.normalized.tsv

echo DoNe
