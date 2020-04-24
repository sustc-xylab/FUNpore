options(echo=F) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

library(plyr)
library(data.table)

krak<-fread("KRAKEN/combined.krak",header=F)
krak<-krak[,2:3]
colnames(krak)<-c("contig","tax_id")

#lineage<-fread("database/lineages-2019-02-20.csv")
lineage<-fread(args[1])
colnames(lineage)[2]<-"kingdom"
lineage<-lineage[,c("tax_id","kingdom","phylum","class","order","family","genus","species")]

kraken<-merge(krak,lineage, by="tax_id")
kraken<-kraken[,-1]


taxator<-fread("classify_bin/contig_taxonomy.tab.modified",fill=T)
rank<-c("kingdom","phylum","class","order","family","genus","species")
if(ncol(taxator)<8){rank<-rank[1:(ncol(taxator)-1)]} # subset rank in case none of the ont read is assigned to species
colnames(taxator)<-c("contig",rank)

# use which kraken as basis to merge in taxator,
# merge in taxa from taxator only when taxator's above level taxa is equal to karen

k<-data.frame(kraken)
t<-data.frame(taxator)
m<-merge(k,t, by="contig", all=T)

for(i in 1:length(rank)){
	x<-paste(rank[i],"x",sep=".") # corespond to k
	y<-paste(rank[i],"y",sep=".") # corespond to t
	lookat<-which(m[,x]=="")
	if(i==1) {m[,x][lookat]<-m[,y][lookat]} 
	else {
		x2<-paste(rank[i-1],"x",sep=".") # corespond to k
		y2<-paste(rank[i-1],"y",sep=".") # corespond to t
		lookat2<-which(m[,x2]==m[,y2])
		lookat3<-intersect(lookat,lookat2)
		if(length(lookat3)>0) {m[,x][lookat3]<-m[,y][lookat3]}
	}
}
m<-m[,1:8]
colnames(m)<-c("contig",rank)
write.table(m,file=args[2],row.name=F,col.name=T, quote=F, sep="\t")

	
