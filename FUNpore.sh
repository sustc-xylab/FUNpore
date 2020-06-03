#!/bin/bash
## FUNpore is designed to do frame-shift correction of nanopore 1D/2D reads (fasta format) 
## and then carried out taxonomic annotation of the corrected reads with taxator-tk and KRAKEN
## Author: Yu XIA - 2020-04-20
## Email: shuixia100@gmail.com
##version 1.0
set -e


#### usage info ####
show_help() {
cat << EOF
Usage: ${0##*/} 
version 1.0
written by: Yu XIA <shuixia100@gmail.com>

arguments:
	-h	display this help 

	-f 	1D.fasta generated by nanopore sequncing as input for FUNpore

	-t 	number of threads used for parallel, default t=1
	
	-c	number of query per parallel run default c=200000, set c=5 or smaller if your average read length is larger than 1M

output files are in directory of FUNpore_nowtime:
	input_framecorrect.fa 	frame-shift corrected nanopore reads
	input_framecorrect.summary.csv	summary of the frame-shift correction 
	input_framecorrect.fa_taxa.tab	phylogenetic assignment of frame-shift corrected reads by taxator-tk(nt based homology search) and KRAKEN (kmer search)


Example usage: 
	step1 create a working directory for FUNpore run and copy the target fasta to the created directory
		mkdir testdir 
		cp test.fa testdir 
	step2 run FUNpore in the directory created
		cd testdir
		bash FUNpore.sh -f ./test.fa -t 20 -c 200000 > FUNpore.log

NOTICE: 
	Make sure there is no previous intermediate files/directory exist in current working directory 
	One given working directory can only be used by one instance of FUNpore run
	use bash not sh
EOF
}

		
####################
# define arguments
####################
# the DIR of FUNpore scirpt
SCRIPT=`realpath $0`
DIR=`dirname $SCRIPT`
nowt=`date +%Y-%m-%d.%H:%M:%S`;

OPTIND=1  # Reset in case getopts has been used previously in the shell.

# set default value:
N_threads="1"
Input_fa=""
cut="200000"


while getopts "t:f:l:s:c:o:h" opt; do
	case "$opt" in
		h|--help)
			show_help
			exit 0
			;;
		t)
			N_threads=$OPTARG
			;;
		f)
			Input_fa=$OPTARG
			;;
		c)
			cut=$OPTARG
			# number of query per parallel run 
			;;
		\?)
			echo "Invalid option: -$OPTARG" >&2
			exit 0
			;;
		'?')
			show_help >&2
			exit 1
			;;
		-?*)
			print 'Warning: Unknown option (ignored) : %s\n' "$1" >&2
			exit 0
			;;
		:)
			echo "Option -$OPTARG requires an argument." >&2
			exit 1
			;;
		*) # default case: if no more options then break out of the loop
			break
			
	esac
done

if [ -z "$Input_fa" ]
then
	echo "No input fasta, -f must be specified"
	exit
fi

shift "$((OPTIND-1))"


export LD_LIBRARY_PATH=${DIR}/bin/lib:$LD_LIBRARY_PATH
export LDFLAGS="-L/${DIR}/bin/lib"


echo "FUNpore is runing using parameters:
------------------------------------------
Input contigs: $Input_fa
Number of threads: $N_threads
Number of query per parallel run: $cut
------------------------------------------
"



# subset the name of the $Input_fa, remove blank
myarray=(`echo $Input_fa| tr "/" " "`)
Input_fa2=${myarray[-1]}

echo "
---------------------------------------------------------------------------------
remove reads shorter than 1kb and rename reads to avoid illegal format in header"
$DIR/bin/fasta.remove.shorter.pl 1000 $Input_fa > ${Input_fa}.1kb+
grep ">" ${Input_fa}.1kb+ | sed 's/>//'> tmp.name1
awk '/^>/{print ">ont" ++i; next}{print}' < ${Input_fa}.1kb+ > tmp
mv tmp ${Input_fa}.1kb+
grep ">" ${Input_fa}.1kb+ | sed 's/>//'> tmp.name2
paste tmp.name1 tmp.name2 > ${Input_fa}.1kb+_original.name
rm -f tmp.name1 tmp.name2

Input_fa2=${Input_fa}.1kb+
 
###############################################
# contig Last against refseq_protein for frame-shift correction
###############################################
echo "
Start Frame-shift correction of $Input_fa2 @ `date +"%Y-%m-%d %T"`
"
out="LAST"
DB_last=`grep "DB_last" ${DIR}/FUNpore_CONFIG | head -1 | sed 's/DB_last=//'`
DB_name=`grep "DB_name" ${DIR}/FUNpore_CONFIG | head -1 | sed 's/DB_name=//'`
# lastdb -Q 0 -P 50 -p -cR01 -v refseq_protein20190808_lastdb ../refseq_protein20190808.fa
if [ ! -d $out ]; then
	mkdir $out;
else
	echo "Warning: $out already exists and will be overwrite."
fi
echo "
LAST alignning agasint refseq_protein database. This step may takes some time please be patient
"
${DIR}/bin/last-983/scripts/parallel-fasta "${DIR}/bin/last-983/src/lastal -f MAF -F 15 -pPAM30 -d 100 -k 2 ${DB_last}/${DB_name} -P $N_threads | ${DIR}/bin/last-983/scripts/maf-convert blast" < ${Input_fa2} > ${out}/${Input_fa2}_last.maf.blast

# The PAM30 scoring scheme finds strong protein similarities:
# -d: Extend a gapless alignment from each initial match, and keep those with score ≥ d. 
# -k: Look for initial matches starting only at every STEP-th position in each query (positions 0, STEP, 2×STEP, etc). This makes lastal faster but less sensitive.

echo "converting last output for frame-shift correction"
grep "Query" ${out}/${Input_fa2}_last.maf.blast | tr -s ' ' '\t' | sed 's/\\/?/g' >  ${out}/${Input_fa2}_last.maf.blast_nanopore.align 

rm -f ${out}/${Input_fa2}_last.maf.blast # taking too much disk space

# seperate every Query=
echo "reordering last output to avoid duplicated entries in results @ `date +"%Y-%m-%d %T"`"
echo "	seperating into single entry by Query= @ `date +"%Y-%m-%d %T"`"

# parallel runing awk 
cd ${out}
# every 20000 Query= cut into subfiles
cat ./${Input_fa2}_last.maf.blast_nanopore.align | parallel --pipe -L20000 -N1 --recstart Query= 'cat >cut_{#}'
# further awk in parallel to cut every Query= into a single file
rm -f tmp.awk.jobs
find ./ -name "cut_*" > tmp.cut
cat tmp.cut| while read line
do
	awk -v line="$line" 'BEGIN{i=0}/Query=\tont*/{i++}{print > "./"line"_split_"i"_nanopore.align"}' $line &
	PID=$!
	echo $PID >> tmp.awk.jobs
done

bash $DIR/bin/monitor.bgpid.sh tmp.awk.jobs 

rm -f tmp.awk.jobs
cat tmp.cut | xargs rm -f
rm -f tmp.cut

cd ..

echo "	further seperating into single entry by ont @ `date +"%Y-%m-%d %T"`" 
rm -rf merge.nanopore.align
mkdir -p merge.nanopore.align
# here cannot simultaneous writing cause simultaneous writing to the same file could result in malformating 
find ${out} -name "cut_*_split_*_nanopore.align" > tmp.align
split -a 4 -d -l 20000 tmp.align tmp.align.cut 
rm -f tmp.combine.jobs
rm -f merge.nanopore.align_ontname

ls tmp.align.cut* | while read line3
do 
  bash $DIR/bin/combine.line.sh $line3 &
  PID=$!
	echo $PID >> tmp.combine.jobs
done

bash $DIR/bin/monitor.bgpid.sh tmp.combine.jobs 
# don't run below when testing
rm -f tmp.combine.jobs
rm -f tmp.align
rm -f tmp.align.cut*
find ${out} -name "cut*" -delete  & # remove the splited files

## combine each query into one ------
echo "	combine ont entries into one and kept only the first Query header in the combined file @ `date +"%Y-%m-%d %T"`"
sort -u ./merge.nanopore.align_ontname > merge.nanopore.align_ontname.uniq
split -a 4 -d  -l 2000 merge.nanopore.align_ontname.uniq merge.nanopore.align_ontname.uniq.cut
rm -rf merge.nanopore.align2
mkdir -p merge.nanopore.align2

# further combine the same ont in parallel 
rm -f tmp.merge.jobs
find merge.nanopore.align -name "*" | cut -d / -f 2 | sort | sed '1d' > tmp.lst
cut -d _ -f 1 tmp.lst > tmp2
find ./ -name "merge.nanopore.align_ontname.uniq.cut*" | while read line
do
  bash $DIR/bin/merge.line.sh $line &
  PID=$!
  echo $PID >> tmp.merge.jobs
done

bash $DIR/bin/monitor.bgpid.sh tmp.merge.jobs
# don't run below commands when testing
rm -f tmp.merge.jobs
rm -f merge.nanopore.align_ontname
rm -f merge.nanopore.align_ontname.uniq
rm -f merge.nanopore.align_ontname.uniq.cut*
rm -f tmp.lst
rm -f tmp2
# fast delete big directory
cd merge.nanopore.align
perl -e 'for(<*>){((stat)[9]<(unlink))}'
cd ..

rm -rf merge.nanopore.align
mv merge.nanopore.align2 merge.nanopore.align

echo "	combine into trunks for subsequent frame-shift correction @ `date +"%Y-%m-%d %T"`"
find ./ -name "tmp.ont*" | xargs rm -f 

ls ./merge.nanopore.align > tmp.ont
split -a 4 -d -l $cut tmp.ont tmp.ont.cut

find ${out} -name "split_*_nanopore.align.modify" | xargs rm -f  # remove intermediates from previous run 
rm -f tmp.merge2.jobs
ls tmp.ont.cut* | while read line
do
	bash $DIR/bin/merge.line2.sh $line $out & 
	PID=$!
	echo $PID >> tmp.merge2.jobs
done

bash $DIR/bin/monitor.bgpid.sh tmp.merge2.jobs
find ./ -name "tmp.ont*" | xargs rm -f  # remove tmp files
rm -f tmp.merge2.jobs

#parallel run samtools in background
echo "subsetting $Input_fa2 according to query @  `date +"%Y-%m-%d %T"`"
samtools faidx $Input_fa2
ls ${out}/split_*_nanopore.align.modify | while read line
do	
	fgrep "Query=" $line | cut -f 2 | sort -u > ${line}_header
	xargs samtools faidx $Input_fa2 < ${line}_header > ${line}_sub.fa &
	PID=$!
	echo $PID >> tmp.samtools.jobs
done

bash $DIR/bin/monitor.bgpid.sh tmp.samtools.jobs
rm -f tmp.samtools.jobs

echo "frame-shift correction in R @ `date +"%Y-%m-%d %T"`"
Rscript ${DIR}/bin/frame_shift.correction.R ${out} $N_threads

echo "combine the splited results into the final result @ `date +"%Y-%m-%d %T"`"
cat ${out}/split_*_nanopore.align.modify_framecorrect.fa > ${Input_fa2}_framecorrect.fa
cat ${out}/split_*_nanopore.align.modify_framecorrect.summary.csv > ${Input_fa2}_framecorrect.summary.csv

echo "remove intermediate files @ `date +"%Y-%m-%d %T"`"
# fast delete big directory
cd merge.nanopore.align
perl -e 'for(<*>){((stat)[9]<(unlink))}'
cd ..
rm -rf merge.nanopore.align # remove intermediate files
find $out -name "*_header" -delete  &
find $out -name "*_sub.fa" -delete  &
find $out -name "*_sub.fa.fai" -delete & 

echo "Done Frame-shift correction @ `date +"%Y-%m-%d %T"`"

###############################################################
# taxonomy annotation of combined.fa by KRAKEN amd taxator-tk
###############################################################
# taxonomy assignment of contig by taxator-tk ########################
# taxator-tk is homology based on blastn search agaisnt nt ####
out="classify_bin"
BLASTDB=`grep "BLASTDB" ${DIR}/FUNpore_CONFIG | head -1 | sed 's/BLASTDB=//'`
TAXDUMP=`grep "TAXDUMP" ${DIR}/FUNpore_CONFIG | head -1 | sed 's/TAXDUMP=//' | sed 's/"//g'`
TAXDUMP2=`echo $TAXDUMP | sed 's/\/nodes.dmp//'`
Query="${Input_fa2}_framecorrect.fa"

echo "
--------------------------------------------------------------------------
Start taxonomy annotation of $Query @ `date +"%Y-%m-%d %T"`
"

echo "Start taxator-tk @ `date +"%Y-%m-%d %T"`"

export TAXATORTK_TAXONOMY_NCBI=$TAXDUMP2

if [ ! -d $out ]; then
        mkdir $out;
else
        echo "Warning: $out already exists and will be overwrite"
fi


echo "
aligning $Query to ${BLASTDB} database with MEGABLAST. please be patient. You may check the progress in  ${out}/megablast_out.raw.tab
"
${DIR}/bin/blastn -task megablast -num_threads $N_threads\
 -db ${BLASTDB}\
 -outfmt '6 qseqid qstart qend qlen sseqid staxids sstart send bitscore evalue nident length'\
 -query ${Query} > ${out}/megablast_out.raw.tab

echo "removing unnecessary lines that lead to bad tax IDs, those without a proper rank"

python ${DIR}/bin/prune_blast_hits.py ${TAXDUMP} ${out}/megablast_out.raw.tab > ${out}/megablast_out.pruned.tab

cat ${out}/megablast_out.pruned.tab | cut -f1,2,3,4,5,7,8,9,10,11,12 > ${out}/megablast_out.tab

echo "making mapping file"
cat ${out}/megablast_out.pruned.tab | cut -f 5,6 > ${out}/mapping.tax


echo "pulling out classifications with taxator using megan-lca algorithm"
cat ${out}/megablast_out.tab | ${DIR}/bin/taxator -a megan-lca -t 0.3 -e 0.00001 -g ${out}/mapping.tax -p 40 > ${out}/predictions.gff3

echo "binning and consolidating classifications for each contig"
sort -k1,1 ${out}/predictions.gff3 | ${DIR}/bin/binner -n classification -i genus:0.6 > ${out}/binned_predictions.txt
mv binning.log ${out}

echo "pulling out full taxonomy path with taxknife"
cat ${out}/binned_predictions.txt | ${DIR}/bin/taxknife -f 2 --mode annotate -s path | grep -v "Could not" | cut -f1,2 > ${out}/contig_taxonomy.tab

# modify for merging with KRAKEN
cat ${out}/contig_taxonomy.tab| sed 's/;/\t/g' > ${out}/contig_taxonomy.tab.modified

echo "
Finish taxator-tk @ `date +"%Y-%m-%d %T"`"


##### taxonomy annotation based on KRAKEN  #################
echo "
Start KRAKEN @ `date +"%Y-%m-%d %T"`"

out="KRAKEN"
KRAKEN_DB=`grep "KRAKEN_DB" ${DIR}/FUNpore_CONFIG | head -1 | sed 's/KRAKEN_DB=//' | sed 's/"//g'`
Query="${Input_fa2}_framecorrect.fa"

if [ ! -d $out ]; then
	mkdir $out;
else
	echo "Warning: $out already exists and will be overwrite"

fi
echo "KRAKEN annotating"
${DIR}/bin/libexec/kraken --db ${KRAKEN_DB} --fasta-input --threads $N_threads --output ${out}/combined.krak $Query

echo "
Finish KRAKEN @ `date +"%Y-%m-%d %T"`"

## merge taxonomy annotation by KRAKEN and taxator-tk ########
echo "merging KRAKEN and taxator-tk annoation in R"
Rscript ${DIR}/bin/combine.kraken.taxator.R ${DIR}/database/lineages-2019-02-20.csv ${Query}_taxa.tab

echo "
Finish taxonomy annotation of $Query @ `date +"%Y-%m-%d %T"`"


echo "
-----------------------------------------------------------------
Saving FUNpore results 
"
out=`echo "${Input_fa}_FUNpore_${nowt}"`
echo "moving results to $out"
if [ ! -d $out ]; then 
	mkdir $out;
else 
	rm -rf $out
	mkdir -f $out

fi

mv ${Input_fa2}* ${out}
mv LAST ${out}
mv classify_bin ${out}
mv KRAKEN ${out}

echo "
done Funpore @ `date +"%Y-%m-%d %T"`
--------------------------------------------------------------------
"
