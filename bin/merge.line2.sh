#!/ubr/bin/bash
cat $1 | while read line2
do
	cat ./merge.nanopore.align/$line2 >> ${2}/split_${1}_nanopore.align.modify
done
