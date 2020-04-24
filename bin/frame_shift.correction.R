
# lastal -f MAF -F 15 ~/db/refseq_protein20170723_lastdb/refseq_protein20170723_lastdb test.fa -P 40 -v > test.fa_refseq_last.maf
# maf-convert blast test.fa_refseq_last.maf > test.fa_refseq_last.maf.blast
# grep "Query" test.fa_refseq_last.maf.blast | tr -s ' ' '\t' | sed 's/\\/?/g' >  test.fa_refseq_last.maf.blast_nanopore.align
library(foreach)
library(doParallel)
library(seqinr)
library(plyr)
library(data.table)
library(Rsamtools)

options(echo=F) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

# testing purpose
# args<-c("LAST", "10")

files<-system(paste("ls ",args[1],"/split_*_nanopore.align.modify",sep=""),intern=T)
fasta.files<-paste(files,"_sub.fa",sep="")
n_threads<-as.numeric(args[2])

cl<-makeCluster(n_threads,outfile="")
registerDoParallel(cl)
foreach(k=1:length(files),.packages=c("plyr","data.table","seqinr","Rsamtools"))%dopar%{
# for(k in 1:length(files)){
	
  # read in nanopore.align and corresponding fasta ----	
	nowtime<-Sys.time()
	cat("start frame-shift correction for ",fasta.files[k], as.character(nowtime),"\n")
		
	fread(files[k],header=F,fill=T,select=1:4)-> d
	d<-as.data.frame(d)
	cat("done read in nanopore.align for", files[k],"\n")
		
	cat("reading in sub.fasta for",fasta.files[k],"\n")
	indexFa(fasta.files[k])   # create an index of file 'foo.fasta'
	fa = FaFile(fasta.files[k])  # reference the fasta file and it's index
	gr = as(seqinfo(fa), "GRanges")
	# fa<-read.fasta(file=fasta.files[k])
	cat("done read in sub.fasta for", fasta.files[k],"\n")
	
	# subset into list ------
	fa.adj<-list()

	cat("subseting into list according to ”Query=“, then subset each Query into item of a list\n") 
	lookat<-grep("Query=",d$V1)
	tmpn<-d[lookat,]$V2

	d.lst<-list()
	if(length(lookat)>1) {
		for(i in 1:length(lookat)){
		if(length(lookat)==1 | i==length(lookat)){
		  tmpr<-d[(lookat[i]+1):nrow(d),]
		}
		else {
		  tmpr<-d[(lookat[i]+1):(lookat[i+1]-1),]
		  tmpr<-tmpr[grep("Query=",tmpr$V1,invert=T),] # to remove the case that two empty query are consecutivly appear and cause format error in subsequent analysis
		}
		d.lst[[i]]<-tmpr
		#if(i%%1000==0 | i == length(lookat)){cat("finish subseting",round(i/length(lookat)*100),"%\n")}
		}
	  names(d.lst)<-tmpn

	  cat("merging d.lst with the same name into one item\n")
	  tmpn.unique<-unique(tmpn)
	  d.lst2<-list()
	  for(c in tmpn.unique){
  		lookat<-which(tmpn==c)
  		tmpd<-rbind.fill(d.lst[lookat])
  		# tmpd<-d.lst[[lookat[1]]]
  		# for(j in 2:length(lookat)){tmpd<-rbind.fill(tmpd,d.lst[[lookat[j]]])}
  		d.lst2[[c]]<-tmpd
  	  }
	  } else {
	  d.lst[[1]]<-d[(lookat[1]+1):nrow(d),]
	  names(d.lst)<-tmpn
	  d.lst2<-d.lst
	  }
	rm(d.lst,d) # to save RAM
	cat("done subset into list according to Query",files[k],as.character(Sys.time()),"\n")
	#save(d.lst2, file=paste(files[k],"_d.lst2.RData",sep=""))
# 	to check if fasta is corresponding with d.lst2
# 	tmp<-sapply(d.lst2,function(x){x<-x[,2:4]
# 	y<-max(as.numeric(c(x[,1],x[,3])))
# 	return(y)})
# 
#   tmp2<-vector()
#   for(i in 1:length(tmp)){tmp2[i]<-length(getSeq(fa, gr[names(tmp)[i]])[[1]])}
# 
#   length(which(tmp<=tmp2))
  
  
	# loop for shift correction ----
	cat("start loop for shift-correction\n")
	result.lst<-list()
	good.f<-list()
	for(g in 1:length(d.lst2)){
		x<-d.lst2[[g]]
		x<-x[,2:ncol(x)]
		colnames(x)<-c("start","alignment","end")
		x$start<-as.numeric(x$start)
		x$end<-as.numeric(x$end)
		x<-as.data.table(x)
		x<-unique(x,by=c("start","end"))
		# x<-as.data.frame(x)
		
		# subset each section of alignment into items of list and reformat 
		lookat<-vector()
		if(nrow(x)>2){
		  for(i in 1:(nrow(x)-1)){
			if((x$start[i+1]!= x$end[i]-1) & (x$start[i+1]!=x$end[i]+1)){lookat[i]<-i} 
			# cat("finish",round(i/nrow(x)*100),"%\n")
		  }
		  i=nrow(x)
		  if((x$start[i]!= x$end[i-1]-1) & x$start[i]!=x$end[i-1]+1 ){lookat[i]<-i} # 正向，反向alignment都处理了
		}
		
		lookat<-lookat[!is.na(lookat)]
		if(length(lookat)==0){lookat<-nrow(x)} else 
		  if(lookat[length(lookat)]!=nrow(x)){lookat<-c(lookat,nrow(x))} # make sure the last item in lookat is the the last line of x 
		
		tmpdfx<-data.frame(start=c(1,lookat+1)[1:length(lookat)],end=lookat)
		x.sub<-apply(tmpdfx,1,function(e){
			y<-x[e[[1]]:e[[2]],]
			m<-paste(y$alignment,collapse="")
			n<-data.frame(start=y$start[1],end=y$end[nrow(y)],alignment=m)
			return(n)
		})
		
		df<-rbind.fill(x.sub)
		
		rm(x.sub) # free up memory
		gc()
		# keep those without frameshift and if such alignment is longer than 0.7 of before correction read length, then this read will be write.out to good.read 
		
		if(length(which(is.na(df$alignment)))>0){
		  df2<-df[grep("\\?|\\/",df$alignment,invert=T),]
		  f<-getSeq(fa, gr[names(d.lst2)[g]])[[1]]
		  # f<-fa[[names(d.lst2)[g]]]
		  if(max(abs(df2$start-df2$end))/length(f)>=0.7){good.f[[names(d.lst2)[g]]]<-f}
		}
		
		df<-df[grep("\\?|\\/",df$alignment),] # remove those alignments don't have frame-shift, before framecorrection
		# identify each position in query that need to do frame-shift correction and then carry out shift-correction -------
		if(nrow(df)>0){
		  df.lst<-list()
		  for(i in 1:nrow(df)){
			e<-df[i,]
		  # df.lst<-apply(df,1,function(e){
			y<-e[[3]]
			y.start<-e[[1]]
			y.end<-e[[2]]
			y<-strsplit(as.character(y),"")[[1]]
			
			if(y.start>y.end){
			  y<-rev(y)
			  loc.1<-grep("\\?|\\/",y)
			  
			  position.lst<-list()
			  for(j in 1:length(loc.1)){
			  	# z<-y[1:(loc.1[j]-1)]
				if((loc.1[j]-1)<1){z<-y[1]} else {z<-y[1:(loc.1[j]-1)]}
				n1=length(grep("\\?",z))
				n2=length(grep("\\/",z))
				n3=length(grep("\\-",z))
				p.start<-y.end+(loc.1[j]-1)*3-3*n3-4*n2-2*n1
				if(p.start>0){position.lst[[j]]<-data.frame(y[loc.1[j]],p.start)} else{position.lst[[j]]<-data.frame(y[loc.1[j]],p.start=1)}
			  }
			} else {
			  loc.1<-grep("\\?|\\/",y)
			  position.lst<-list()
			  for(j in 1:length(loc.1)){
				# z<-y[1:(loc.1[j]-1)]
				if((loc.1[j]-1)<1){z<-y[1]} else {z<-y[1:(loc.1[j]-1)]}
				n1=length(grep("\\?",z))
				n2=length(grep("\\/",z))
				n3=length(grep("\\-",z))
				p.start<-y.start+(loc.1[j]-1)*3-3*n3-4*n2-2*n1
				if(p.start>0){position.lst[[j]]<-data.frame(y[loc.1[j]],p.start)} else{position.lst[[j]]<-data.frame(y[loc.1[j]],p.start=1)}
			  }
			}
			p.df<-rbind.fill(position.lst)
			
			# return(p.df)
			df.lst[[i]]<-p.df
			if(nrow(df)>100000 & (i==1|i%%10000==0)){cat("finish searching df",round(i/nrow(df)*100),"%\n")}
			}
		  
		  tmp<-arrange(rbind.fill(df.lst),p.start)
		  colnames(tmp)<-c("shift.type","position")
		  
		  # frame-shift tha appear in multiple alignments, those frame-shifts are considered as confirmed
		  lookat<-duplicated(tmp)
		  if(length(which(lookat))==0){confirmed.p<-data.frame(shift.type=NA,position=NA,status=NA)} else{lookat<-duplicated(tmp)
		  confirmed.p<-tmp[duplicated(tmp),]
		  confirmed.p<-confirmed.p[!duplicated(confirmed.p),]
		  confirmed.p<-confirmed.p[!duplicated(confirmed.p$position),]
		  confirmed.p$status<-"confirmed"}
		  
		  # frame-shift that only report in one alignments
		  candidate.p<-tmp[which(!tmp$position %in% confirmed.p$position),]
		  if(nrow(candidate.p)==0){candidate.p<-data.frame(shift.type=NA,position=NA,status=NA)}else{candidate.p$status<-"candidatus"}
		  
		  # all frame-shift
		  all.p<-rbind(confirmed.p,candidate.p)
		  all.p<-all.p[which(!is.na(all.p$shift.type)),]
		  all.p<-arrange(all.p,position)
		  
		  result.p<-list(df,all.p)
		  
		  # carry out frame-shift correction of each questy ------------------------
		  # f<-fa[[names(d.lst2)[g]]]
		  f<-getSeq(fa, gr[names(d.lst2)[g]])[[1]]
		  # names(f)<-names(d.lst2)[g]
		  
		  # loop to do the frame-shift correction, only do frame-shift correction on those confirmed frame-shifts
		  all.p.adj<-confirmed.p
		  f.adj<-f
		  if(length(which(is.na(all.p.adj$shift.type)))!=nrow(all.p.adj)){
			for(i in 1:nrow(all.p.adj)){
			  if(all.p.adj$shift.type[i]=="/"){
				# nt at this position should be repeat once
				tmpp<-all.p.adj$position[i]-1
				
				# repeat this nt in f.adj
				f.adj<-c(f.adj[1:tmpp],f.adj[tmpp],f.adj[(tmpp+1):length(f.adj)])
				# increase one nt for all the position since one nt is added in the sequence
				all.p.adj$position<-all.p.adj$position+1
			  } else {
				# nt at this position shoulbe be deleted from the sequence
				tmpp<-all.p.adj$position[i]
				
				# delete this nt in f.adj
				if(tmpp>1){f.adj<-c(f.adj[1:(tmpp-1)],f.adj[(tmpp+1):length(f.adj)])} else {f.adj<-c(f.adj[(tmpp+1):length(f.adj)])}
				# decrease one nt for all the position sicne one nt is deleted from the sequence 
				all.p.adj$position<-all.p.adj$position-1
			  }
			  #cat("finish correction",round(i/nrow(all.p.adj)*100),"%\n")
			}
			
			# result.fa<-lapply(f.adj,toupper)
			result.fa<-f.adj
			}else{result.fa<-NULL}
		  result<-list(result.p,result.fa)
		} else {result<-list(NULL,NULL)}

		result.lst[[g]]<-result
		if((length(d.lst2)>=1000 & (g==1 | g%%50==0 | g == length(d.lst2)))| (g==1 | g%%2==0 | g == length(d.lst2))) {cat("finish correction of",files[k],round(g/length(d.lst2)*100),"%\n")}
	}
	cat("done loop for shift-correction for",files[k],as.character(Sys.time()),"\n")
	#save(result,file=paste(files[k],"_result.RData",sep=""))
	
	# write out frame-shift corrected fasta---------------------------------------------------
	fa.adj<-lapply(result.lst,"[[",2)
	names(fa.adj)<-names(d.lst2)
	lookat<-which(sapply(fa.adj,function(x) !is.null(x)))
	cat("writing out the frame-shift corrected fasta of ",files[k],"\n")
	write.fasta(fa.adj[lookat],
				names = names(d.lst2)[lookat],
				file.out=paste(files[k],"_framecorrect.fa",sep=""),
				open="w")
	cat("done write out the frame-shift corrected fasta of ",files[k], as.character(Sys.time()),"\n")
	
	if(length(good.f)==0){cat("there is no good fasta for", files[k],"\n")} else {
	  lookat2<-which(sapply(good.f,function(x) !is.null(x)))
	  cat("writing out the frame-shift corrected fasta of ",files[k],"\n")
	  write.fasta(good.f[lookat2],
	              names = names(good.f[lookat2]),
	              file.out=paste(files[k],"_good.fa",sep=""),
	              open="w")
	  cat("done write out the good fasta of ",files[k], as.character(Sys.time()),"\n")
	}
  
	
	# wirte out frame-shift correction summary ------
	cat("writing out summary file for", files[k]," frame-shift correction\n")
	d.lstf<-lapply(result.lst,"[[",1)
	df.lst<-lapply(d.lstf,"[[",1)
	p.lst<-lapply(d.lstf,"[[",2)

	correction.summary<-matrix(nrow=0,ncol=6)
	for(g in 1:length(d.lst2)){
		x<-p.lst[[g]]
		tmpq<-names(d.lst2)[g]
		f<-getSeq(fa, gr[tmpq])[[1]]
		# # fa<-read.fasta(file=fasta.files[k])
		# f<-fa[[names(d.lst2)[g]]]
		b.fa.len=length(f)
		a.fa.len=length(fa.adj[[tmpq]])
		if(is.null(x)){n.correction=0
		n.bs=0
		n.fs=0} else{
		  t<-aggregate(position~status+shift.type,x,length)
		  t<-t[t$status=="confirmed",]
		  if(nrow(t)>0){n.correction=sum(t$position)}else{n.correction=0}
		  if(length(which(t$shift.type=="?"))>0){n.bs= t$position[which(t$shift.type=="?")]}else{n.bs=0}
		  if(length(which(t$shift.type=="/"))>0){n.fs= t$position[which(t$shift.type=="/")]}else{n.fs=0}
		}
		
		y<-t(as.matrix(c(tmpq,n.correction,n.bs,n.fs,b.fa.len,a.fa.len),bycol=T))
		correction.summary<-rbind(correction.summary,y)
		#if(g%%100==0 | g == length(d.lst2)){cat("finish ",round(g/length(d.lst2)*100),"%\n")}
	}

	colnames(correction.summary)<-c("query","n.correction","n.bs","n.fs","before.correction.len","after.correction.len")
	write.csv(correction.summary,file=paste(files[k],"_framecorrect.summary.csv",sep=""),row.names=F)
	
	correction.summary<-data.frame(correction.summary)
	correction.summary[,2:6]<-apply(correction.summary[,2:6],2, function(x) as.numeric(as.character(x)))
	n<-sum(correction.summary$n.correction)/sum(correction.summary$before.correction.len)*1000
	cat("there are",n,"frame-shifts corrected per kbp of ",files[k],"\n\n")
	cat("done frame-shift correction for ",fasta.files[k],as.character(Sys.time()),"\n")
	
	gc() # to free up system RAM
}


stopCluster(cl)

























