#!/bin/bash
# script to combine cutted align according to query
i=1
max=`wc -l $1 | cut -f 1`

cat $1 | while read line
do
  line2=`echo $line | cut -d "/" -f 2`
  name=`head -1 $line | sed 's/Query=\t//'`
  cat $line > ./merge.nanopore.align/${name}_${line2}
  echo "$name" >> ./merge.nanopore.align_ontname
  
  # # scripts for progress bar
  # echo -ne "\r Finished $i/$max $1"
  # echo $((i+=1)) >/dev/null
  
done
echo "finish $1"
	

