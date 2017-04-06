#!/bin/bash

vcftools \
	--gzvcf Partial_Complete_F2.vcf.gz \
	--out F2_filtered \
	--positions parent_postions.tab \
	--remove-indels \
	--min-alleles 2 \
	--max-alleles 2 \
	--minQ 40 \
	--minDP 10 \
	--maxDP 1000 \
	--recode \
	--recode-INFO-all

mv F2_filtered.recode.vcf F2_filtered.vcf
gzip F2_filtered.vcf
