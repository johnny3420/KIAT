#!/bin/bash

files=`ls -d *.star.dir`
echo $files
echo "Alignment	Unique	Multi	Unmapped" >  mapping.stats
for i in $files
	do
	echo $i | sed 's/\.star\.dir//' > temp.name
	sed -n '10p' ${i}/Log.final.out | awk '{print $6}' | sed 's/%//' > temp.unique
	sed -n '25p' ${i}/Log.final.out | awk '{print $9}' | sed 's/%//' > temp.multi.1
	sed -n '27p' ${i}/Log.final.out | awk '{print $10}' | sed 's/%//' > temp.multi.2
	paste temp.multi.1 temp.multi.2 | awk '{print ($1 + $2)}' > temp.multi.3
	paste temp.multi.3 temp.unique | awk '{print (100 - $1 - $2)}' > temp.unmapped
	paste temp.name temp.unique temp.multi.3 temp.unmapped >> mapping.stats
	rm temp*
	done
