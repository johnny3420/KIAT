#!/bin/bash

cd ~/KIAT/STAR/50bp/single_end/50bp-05-11-Final-8ul_1.fq.star.dir/ &&
grep "^@SN" Unmapped.out.mate1 |cut -f 1 > unMapped.names.mate1 &&
pwd &&
echo done &&
cd ~/KIAT/STAR/50bp/single_end/50bp-05-11-Final-8ul_2.fq.star.dir/ &&
grep "^@SN" Unmapped.out.mate1 |cut -f 1 > unMapped.names.mate1 &&
pwd &&
echo done &&
cd ~/KIAT/STAR/100bp/single_end/100bp-05-11-Final-8ul_1.fq.star.dir/ &&
grep "^@SN" Unmapped.out.mate1 |cut -f 1 > unMapped.names.mate1 &&
pwd &&
echo done &&
cd ~/KIAT/STAR/100bp/single_end/100bp-05-11-Final-8ul_2.fq.star.dir/ &&
grep "^@SN" Unmapped.out.mate1 |cut -f 1 > unMapped.names.mate1 &&
pwd &&
echo done &&
cd ~/KIAT/STAR/100bp/paired_end/100bp-05-11-Final-8ul.star.dir/ &&
grep "^@SN" Unmapped.out.mate1 |cut -f 1 > unMapped.names.mate1 &&
grep "^@SN" Unmapped.out.mate2 |cut -f 1 > unMapped.names.mate2 &&
cat unMapped.names.* > unMapped.names
pwd &&
echo done &&
cd


