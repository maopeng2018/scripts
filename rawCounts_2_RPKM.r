setwd("D:\\Projects\\Basidiomycetes_genomes\\transcriptome_data\\selected_data\\Pycnoporus_coccineus_BRFM310\\GSE74234_RAW");

file_len<- read.delim("Pyccol_exon_length.txt",header=T,sep="\t")

file_count<- read.delim("raw_counts_Pycco1.txt",header=T,sep="\t")
colnames(file_len)<- c("GeneName","Len")

combine <- merge(file_len,file_count,by="GeneName") 

oneB<-as.double(10^9)
RPKMs=combine;

  
for(i in 3:ncol(combine)){
total_count<- sum(combine[,i])
 
for(j in 1:nrow(RPKMs)){
  RPKMs[j,i]<- (oneB*as.double(combine[j,i]))/(as.double(total_count)*combine[j,2]);
 }

}

write.table(RPKMs,file="RPKMs.txt",sep="\t", col.names = T, row.names = F)
write.table(combine,file="combine.txt",sep="\t", col.names = T, row.names = F)

