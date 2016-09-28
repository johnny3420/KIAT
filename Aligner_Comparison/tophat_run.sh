#!/bin/bash

tophat2 -p 4 -G ../Tests/test.Brassica_napus.annotation_v5.gff3 -o 2_100bp.paired Brassica_napus_v4.1.chromosomes ../raw_data/Flowering/100bp_2_1.fq ../raw_data/Flowering/100bp_2_2.fq &&
tophat2 -p 4 -G ../Tests/test.Brassica_napus.annotation_v5.gff3 -o 2_100bp.single Brassica_napus_v4.1.chromosomes ../raw_data/Flowering/100bp_2_1.fq &&
tophat2 -p 4 -G ../Tests/test.Brassica_napus.annotation_v5.gff3 -o 2_50bp.paired Brassica_napus_v4.1.chromosomes ../raw_data/Flowering/50bp_2_1.fq ../raw_data/Flowering/50bp_2_2.fq &&
tophat2 -p 4 -G ../Tests/test.Brassica_napus.annotation_v5.gff3 -o 2_50bp.single Brassica_napus_v4.1.chromosomes ../raw_data/Flowering/50bp_2_1.fq &&
echo DoNe
