#!/bin/bash
# script to combine cutted align according to query
i=1
max=`wc -l $1 | cut -f 1`

cat $1 | while read line
do
  #cat ./merge.nanopore.align/${line}_cut_*_split_*_nanopore.align >> ./merge.nanopore.align2/${line} # too slow
  # find merge.nanopore.align -name "$line" -exec cat '{}' >> ./merge.nanopore.align2/$line \; #even slower 
  start=`fgrep -n -Fw $line tmp2 | cut -d : -f 1 | head -1`
  end=`fgrep -n -Fw $line tmp2 | cut -d : -f 1 | tail -1`
  sed -n "${start},${end}p" tmp.lst | while read line2
  do
	cat merge.nanopore.align/$line2 >> merge.nanopore.align2/$line
	sed  -i '1!{/^Query=\t*/d}' merge.nanopore.align2/$line
  done
done
echo "finish $1"