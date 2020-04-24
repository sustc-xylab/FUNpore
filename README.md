# FUNpore

FUNpore is designed to do frame-shift correction of nanopore 1D/2D reads (fasta format) and  then carried out taxonomic annotation of the corrected reads with taxator-tk and KRAKEN
Author: Yu XIA - 2020-04-20
Email: shuixia100@gmail.com
version 1.0


###### Installation ########
# pre-requisites for FUNpore 
	ruby 2.3.1p112
	samtools 1.9
	R and library: foreach, doParallel, seqinr, plyr, data.table, Rsamtools
	parallel

# download the lastest FUNpore
	
	git clone https://github.com/sustc-xylab/FUNpore.git
	
	cd FUNpore

Once download FUNpore package, all needed analysis is wrapped up in one executable named FUNpore.sh 
The FUNpore_CONFIG contains the PATH for database required for FUNpore, this file should always be stored in the same directory with FUNpore.sh. 
Before runing FUNpore, users should modify FUNpore_config with their specific database PATH

# tools supported out-of-the-box are
	lastal, there is one copy of last-938 included in FUNpore package
	blast+, there is one copy of BLAST 2.2.31+ included in FUNpore package
	taxator-tk, there is one copy of taxator-tk included in FUNpore package
	KRAKEN, there is one copy of KRAKEN 0.10.6-unreleased included in FUNpore package


	
###### Database to prepare before FUNpore run ######
	1. LASTAL database: we recommend run lastal against NCBI refseq_protein database, you may download the blast-formatted refseq_protein database fron NCBI and then extract fasta from the blast-formatted database downloaded, finally build last index with refseq_protein
		1. download the blast-formatted refseq_protein database
		
			perl  update_blastdb.pl --passive --decompress refseq_protein
		
		2. extracted the fasta from the blast-formatted refseq_protein database
		
			blastdbcmd -db refseq_protein -dbtype prot -entry all -outfmt "%f"  -out refseq_protein.fa 
	
		3. build last index with the fasta file extracted
		
			lastdb -Q 0 -P 50 -p -cR01 -v refseq_protein_lastdb refseq_protein.fa
		
		4.specify the path and name of the last database in FUNpore_CONFIG 
		
	2. NCBI nt database for BLAST+ and nodes.dmp of corresponding taxonomy 
		1.download NCBI preformatted nt database with:
			
			perl  update_blastdb.pl --passive --decompress nt
			
		specify the name of the nt database in FUNpore_CONFIG
		
		2. download NCBI taxonomy with：
			wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
		
		Place the files names.dmp and nodes.dmp in a folder and specify its path in FUNpore_CONFIG. keep in mind that the taxonomy files are modified on a regular basis.  
		
	3. default database of KRAKEN
		To create the standard Kraken database, you can use the following command:
		
			FUNpore/bin/libexec/kraken/kraken-build --standard --threads 24 --db $DBNAME
		
		Replace "$DBNAME" with your preferred database name/location.
		
		specify the path of the KRAKEN database in FUNpore_CONFIG
	4. lineage databse converted from NCBI taxonomy
		wget https://gitlab.com/zyxue/ncbitax2lin-lineages/blob/master/lineages-2019-02-20.csv.gz
		gunzip lineages-2019-02-20.csv.gz
		mv lineages-2019-02-20.csv $PATH_to_FUNpore/database

##### Using FUNpore ########
There are two steps to run FUNpore as below. Please NOTICE, since FUNpore will overwrite intermediate files from previous run, as a result, each working directory can only be used for ONE instance of FUNpore run. Additionally, please use bash instead of sh to initiate FUNpore.

	step1 create a working directory for FUNpore run and copy the target fasta to the created directory
		mkdir testdir 
		cp test.fa testdir 
		
	step2 run FUNpore in the directory created, $PATH is the FUNpore installing directory 
		cd testdir
		bash $PATH_to_FUNpore/FUNpore.sh -f ./test.fa -t 20 -c 200000 > FUNpore.log



*Citation:*
if you use FUNpore in your nanopore dataset analysis please cite:
Xia, Yu, An-Dong Li, Yu Deng, Xiao-Tao Jiang, Li-Guan Li, and Tong Zhang. MinION Nanopore Sequencing Enables Correlation between Resistome Phenotype and Genotype of Coliform Bacteria in Municipal Sewage. Frontiers in Microbiology 2017



