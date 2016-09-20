### Setup work space
```
cd
mkdir KIAT
cd KIAT/
mkdir raw_data reference New_STAR Kallisto
```
### Gathering files
```
cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/2016_summer/raw_data
### Using fastq files suggested by Ruijuan
cp 50bp-05-11-Final-8ul_1.fq 50bp-05-11-Final-8ul_2.fq KIAT/raw_data/
cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/ruijuanli/Reference/B.napus
cp star_genome Brassica_napus_v4.1.chromosomes.fa Brassica_napus.annotation_v5.gff3 KIAT/reference/
```

### Trimming fq files from 100bp to 50bp
```
cd KIAT/raw_data
fastx_trimmer -l 50 -i 05-11-Final-8ul_1.fq -o 50bp-05-11-Final-8ul_1.fq
fastx_trimmer -l 50 -i 05-11-Final-8ul_2.fq -o 50bp-05-11-Final-8ul_2.fq
```

###Quick gene count with Kallisto
```
cd KIAT/Kallisto
mkdir paired_100bp	paired_50bp	single_100bp	single_50bp
kallisto quant --plaintext -i ../reference/Brassica_napus.annotation_v5.cds.19.kai -o paired_50bp/ ../raw_data/50bp-05-11-Final-8ul_1.fq ../raw_data/50bp-05-11-Final-8ul_2.fq
kallisto quant --plaintext -i ../reference/Brassica_napus.annotation_v5.cds.19.kai -o paired_100bp/ ../raw_data/05-11-Final-8ul_1.fq ../raw_data/05-11-Final-8ul_2.fq
kallisto quant --single --plaintext -s 50 -l 250 -i ../reference/Brassica_napus.annotation_v5.cds.19.kai -o single_50bp/ ../raw_data/50bp-05-11-Final-8ul_1.fq
kallisto quant --single --plaintext -s 50 -l 250-i ../reference/Brassica_napus.annotation_v5.cds.19.kai -o single_100bp/ ../raw_data/05-11-Final-8ul_1.fq
combine_counts.sh #combines all count files into 2 files
```

###Using R to edit gff3 file
```
library(rtracklayer)
gff <- import.gff3("Brassica_napus.annotation_v5.gff3")
gff
gff$gene_id <- ifelse(is.na(gff$ID),gff$Parent,gff$ID)
export(gff,"Brassica_napus.annotation_v5.gff3",format="gff3")
```

###Running Paired End and Single End STAR on 50bp and 100bp reads
####Star Version 2.5.2b on Whitney, not using Coloma since it's on version 2.4.1d and lacks enough Memory
```
cd KIAT/New_star
mkdir 100bp.paired.star.dir	50bp.paired.star.dir 100bp.single.star.dir 50bp.single.star.dir
cd 100bp.paired.star.dir/
STAR --genomeDir ../../reference/star_genome --sjdbGTFfile ../../Tests/test.Brassica_napus.annotation_v5.gff3 --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfeatureExon CDS --sjdbOverhang 99 --readFilesIn ../../raw_data/05-11-Final-8ul_1.fq ../../raw_data/05-11-Final-8ul_2.fq --outSAMtype BAM Unsorted SortedByCoordinate --outSAMunmapped Within --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic --runThreadN 9
cd ../50bp.paired.star.dir/
STAR --genomeDir ../../reference/star_genome --sjdbGTFfile ../../Tests/test.Brassica_napus.annotation_v5.gff3 --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfeatureExon CDS --sjdbOverhang 49 --readFilesIn ../../raw_data/50bp-05-11-Final-8ul_1.fq ../../raw_data/50bp-05-11-Final-8ul_2.fq --outSAMtype BAM Unsorted SortedByCoordinate --outSAMunmapped Within --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic --runThreadN 9
cd 100bp.single.star.dir/
STAR --genomeDir ../../reference/star_genome --sjdbGTFfile ../../Tests/test.Brassica_napus.annotation_v5.gff3 --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfeatureExon CDS --sjdbOverhang 99 --readFilesIn ../../raw_data/05-11-Final-8ul_1.fq --outSAMtype BAM Unsorted SortedByCoordinate --outSAMunmapped Within --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic --runThreadN 9
cd ../50bp.single.star.dir/
STAR --genomeDir ../../reference/star_genome --sjdbGTFfile ../../Tests/test.Brassica_napus.annotation_v5.gff3 --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfeatureExon CDS --sjdbOverhang 49 --readFilesIn ../../raw_data/50bp-05-11-Final-8ul_1.fq --outSAMtype BAM Unsorted SortedByCoordinate --outSAMunmapped Within --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic --runThreadN 9
```
Attempted on Coloma and got segmentation fault 11. After reading up on fault, may have been caused due to running out of memory on Coloma. Worked on Whitney. 

###Combining star read counts into one file
```
awk '{print $1 "\t" $2}' ReadsPerGene.out.tab| sed 1,4d | sed '1 i\ID 'single_50bp'' > ../Read.counts/single_50bp.counts
awk '{print $1 "\t" $2}' ReadsPerGene.out.tab| sed 1,4d | sed '1 i\ID 'single_100bp'' > ../Read.counts/single_100bp.counts
awk '{print $1 "\t" $2}' ReadsPerGene.out.tab| sed 1,4d | sed '1 i\ID 'paired_50bp'' > ../Read.counts/paired_50bp.counts
awk '{print $1 "\t" $2}' ReadsPerGene.out.tab| sed 1,4d | sed '1 i\ID 'paired_100bp'' > ../Read.counts/paired_100bp.counts
paste *.counts | awk '{ for (i=1;i<=NF;i+=2) {$i=""} print $0}' | sed 's/ /\t/g' > reads
paste *.counts | awk '{print $1}' > ID
paste ID reads > star.read.count.tsv
```

###Planned running cufflinks on output
```
cd 50bp.paired.star.dir/
mkdir cufflink_output      
cufflinks -u -p 9 -o cufflink_output --library-type fr-secondstrand -g ../../Tests/test.Brassica_napus.annotation_v5.gff3 -b ../../reference/Brassica_napus_v4.1.chromosomes.fa Aligned.sortedByCoord.out.bam
cd ../101bp.paired.star.dir/
mkdir cufflink_output
cufflinks -u -p 16 -o cufflink_output --library-type fr-secondstrand -g ../../reference/Brassica_napus.annotation_v5.gtf -b ../../Brassica_napus_v4.1.chromosomes.fa Aligned.sortedByCoord.out.bam
```
