#!/bin/bash

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/De_novo_Assembly/Alignments/Da-Ol-1/Mate_1

Trinity --seqType fq --single 1.Unmapped.out.mate1,2.Unmapped.out.mate1,3.Unmapped.out.mate1,4.Unmapped.out.mate1,All1_Cho_2.Unmapped.out.mate1,All1_Cho_3.Unmapped.out.mate1,All1_Chu_2.Unmapped.out.mate1,All1_Chu_3.Unmapped.out.mate1,All1_Gae_2.Unmapped.out.mate1,All1_Gae_3.Unmapped.out.mate1,All1_Hu_2.Unmapped.out.mate1,All1_Hu_3.Unmapped.out.mate1,All1_Yeong2_160621.Unmapped.out.mate1,All1_Yeong3_160621.Unmapped.out.mate1 --SS_lib_type F --max_memory 50G --CPU 4 --output trinity_out_unmapped_mate1

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/De_novo_Assembly/Alignments/Da-Ol-1/Mate_2

Trinity --seqType fq --single 1.Unmapped.out.mate2,2.Unmapped.out.mate2,3.Unmapped.out.mate2,4.Unmapped.out.mate2,All1_Cho_2.Unmapped.out.mate2,All1_Cho_3.Unmapped.out.mate2,All1_Chu_2.Unmapped.out.mate2,All1_Chu_3.Unmapped.out.mate2,All1_Gae_2.Unmapped.out.mate2,All1_Gae_3.Unmapped.out.mate2,All1_Hu_2.Unmapped.out.mate2,All1_Hu_3.Unmapped.out.mate2,All1_Yeong2_160621.Unmapped.out.mate2,All1_Yeong3_160621.Unmapped.out.mate2 --SS_lib_type R --max_memory 50G --CPU 4 --output trinity_out_unmapped_mate2

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/De_novo_Assembly/Alignments/Da-Ae/Mate_1

Trinity --seqType fq --single 05-11-Final-8ul.Unmapped.out.mate1,5.Unmapped.out.mate1,6.Unmapped.out.mate1,8.Unmapped.out.mate1,Ae_Cho_2.Unmapped.out.mate1,Ae_Cho_3.Unmapped.out.mate1,Ae_Chu_2.Unmapped.out.mate1,Ae_Gae_2.Unmapped.out.mate1,Ae_Gae_3.Unmapped.out.mate1,Ae_Hu_2.Unmapped.out.mate1,Ae_Hu_3.Unmapped.out.mate1,Ae_Yeong_2.Unmapped.out.mate1,Ae_Yeong_3.Unmapped.out.mate1 --SS_lib_type F --max_memory 50G --CPU 4 --output trinity_out_unmapped_mate1

cd /Network/Servers/avalanche.plb.ucdavis.edu/Volumes/Mammoth/Users/johndavis/KIAT/De_novo_Assembly/Alignments/Da-Ae/Mate_2

Trinity --seqType fq --single 05-11-Final-8ul.Unmapped.out.mate2,5.Unmapped.out.mate2,6.Unmapped.out.mate2,8.Unmapped.out.mate2,Ae_Cho_2.Unmapped.out.mate2,Ae_Cho_3.Unmapped.out.mate2,Ae_Chu_2.Unmapped.out.mate2,Ae_Gae_2.Unmapped.out.mate2,Ae_Gae_3.Unmapped.out.mate2,Ae_Hu_2.Unmapped.out.mate2,Ae_Hu_3.Unmapped.out.mate2,Ae_Yeong_2.Unmapped.out.mate2,Ae_Yeong_3.Unmapped.out.mate2 --SS_lib_type R --max_memory 50G --CPU 4 --output trinity_out_unmapped_mate2

echo DOOOOOONNNEEEE
