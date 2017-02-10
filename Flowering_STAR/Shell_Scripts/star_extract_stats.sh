#!/bin/bash

echo Collecing names of all directories which end in .star.dir

files=`ls -d *.star.dir`
#echo $files
echo `wc -w <<< $files` Directories found

awk 'BEGIN { print "Sample\tNumber_Input_Reads\tNumber_Unique_Mapped\tPercent_Unique_Mapped\tNumber_Multi_Mapped\tPercent_Multi_Mapped\tNumber_Too_Many_Multi_Mapped\tPercent_Too_Many_Multi_Mapped\tPercent_Unmapped_Mismatches\tPercent_Unmapped_Too_Short\tPercent_Unmapped_Other\tNumber_Chimeric_Reads\tPercent_Chimeric_Reads" }' > Star_Stats.tab

for i in $files
    do
    echo Extracting mapping stats from Log.final.out file in ${i}
    if [ -f ${i}/Log.final.out ]
        then
            echo ${i} | sed 's/\.star\.dir//' > star.temp.0
            sed -n '6p' ${i}/Log.final.out | awk '{print $6}' > star.temp.1
            sed -n '9p' ${i}/Log.final.out | awk '{print $6}' > star.temp.2
            sed -n '10p' ${i}/Log.final.out | awk '{print $6}' | sed 's/%//' > star.temp.3
            sed -n '24p' ${i}/Log.final.out | awk '{print $9}' > star.temp.4
            sed -n '25p' ${i}/Log.final.out | awk '{print $9}' | sed 's/%//' > star.temp.5
            sed -n '26p' ${i}/Log.final.out | awk '{print $10}' > star.temp.6
            sed -n '27p' ${i}/Log.final.out | awk '{print $10}' | sed 's/%//' > star.temp.7
            sed -n '29p' ${i}/Log.final.out | awk '{print $9}' | sed 's/%//' > star.temp.8
            sed -n '30p' ${i}/Log.final.out | awk '{print $8}' | sed 's/%//' > star.temp.9
            sed -n '31p' ${i}/Log.final.out | awk '{print $7}' | sed 's/%//' > star.temp.10
            sed -n '33p' ${i}/Log.final.out | awk '{print $6}' > star.temp.11
            sed -n '34p' ${i}/Log.final.out | awk '{print $6}' | sed 's/%//' > star.temp.12
            paste star.temp.0 star.temp.1 star.temp.2 star.temp.3 star.temp.4 star.temp.5 star.temp.6 star.temp.7 star.temp.8 star.temp.9 star.temp.10 star.temp.11 star.temp.12 > star.temp
            cat star.temp >> Star_Stats.tab
    else
        echo Error, Log.final.out not found in ${i}
    fi
done

rm star.temp*
echo Star_Stats.tab has been created
