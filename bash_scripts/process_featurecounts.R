
# This script will accept a directory prefix as an argument and 
# (1) plot pseudoalignment quality measurements
# (2) merge the featurecounts abundance results for all samples

CLargs<-commandArgs(TRUE)

dir_name<-CLargs[1]


# Plotting featureCounts quality results

temp<-read.csv(paste(dir_name,'/all_run_info.txt',sep=''),sep='\t',row.names=1)
featurecounts_quality_dat<-as.data.frame(t(temp))

pdf(paste(dir_name,'/featurecounts_results.pdf',sep=''))
barplot(featurecounts_quality_dat$Assigned,names.arg=dimnames(featurecounts_quality_dat)[[1]],main='# of aligned reads annotated',cex.names=0.75,las=2)

barplot(featurecounts_quality_dat$Assigned/rowSums(featurecounts_quality_dat),names.arg=dimnames(featurecounts_quality_dat)[[1]],main='% of aligned reads annotated',cex.names=0.75,las=2,ylim=c(0,1))
dev.off()



# Merging featureCounts count results

fnames<-dir(dir_name,pattern='.txt',recursive=TRUE)
fnames_discard<-c('all_run_info.txt',dir(dir_name,pattern='.summary',recursive=TRUE))
fnames<-fnames[!fnames%in%fnames_discard]
dir_fnames<-paste(dir_name,'/',fnames,sep='')
colnames4dat<-unlist(strsplit(unlist(lapply(strsplit(fnames,split='featurecounts_'),function(x){x[2]})),split='/'))

for(i in 1:length(dir_fnames)){
  dat_i<-read.table(dir_fnames[i],sep='\t',header=TRUE)
  if(i==1){
    length_dat<-data.frame(dat_i[,1], dat_i[,6])
    counts_dat<-data.frame(dat_i[,1], dat_i[,7])
    names(length_dat)<-c('gene',colnames4dat[i])
    names(counts_dat)<-c('gene',colnames4dat[i])
  }
  if(i>1){
    temp_length<-data.frame(dat_i[,1], dat_i[,6])
    temp_counts<-data.frame(dat_i[,1], dat_i[,7])
    names(temp_length)<-c('gene',colnames4dat[i])
    names(temp_counts)<-c('gene',colnames4dat[i])

    length_dat<-merge(length_dat,temp_length,by="gene",all=TRUE)
    counts_dat<-merge(counts_dat,temp_counts,by="gene",all=TRUE)
  }
}

write.table(length_dat,file=paste(dir_name,'/gene_length.tsv',sep=''),row.names=FALSE,quote=FALSE,sep='\t')
write.table(counts_dat,file=paste(dir_name,'/gene_count.tsv',sep=''),row.names=FALSE,quote=FALSE,sep='\t')



