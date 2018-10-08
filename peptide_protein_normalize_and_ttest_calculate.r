setwd("D:\\Projects\\proteomics_EMSL\\Master_sheets\\Extracellular_ttest_uniquePep");

data=read.table("Combined_Extracellular_Master_Sheet_ZeroFill_unique_pep_only_no_decoy.txt",sep="\t",head=T);

colnames(data)

x2=data[,c(3,5:ncol(data))];


proteinRollup= aggregate(. ~ Reference, data=x2, FUN=sum)





x=as.matrix(proteinRollup[,2:ncol(proteinRollup)])   ### be sure which colum start to have the real instensity values, here is colum-5  
adjust=mean(colMeans(x));

x3=x;
for(i in 1:ncol(x)){
x3[,i]=x[,i]*(adjust/mean(x[,i]))
}
x4=log2(x3+1);

name2=colnames(proteinRollup)[2:ncol(proteinRollup)]
name3=paste(name2,"_meanNormalize",sep="")

name4=paste(name2,"_meanNormalize_Log",sep="")
name5=c(colnames(proteinRollup),name3,name4)

proteinRollup2=cbind(proteinRollup,x3,x4)
colnames(proteinRollup2)=name5

write.table(proteinRollup2,"proteinRollup_intensity.txt",sep="\t",col.names = NA)







comparison=read.table("comparisons.txt",sep="\t",head=T);
comps=comparison$comparisons



foldChange <- function (x) {return(tryCatch(sum(2^x[exp]-1)/sum(2^x[cont]-1), error=function(e) NA))}
## foldChange <- function (x) {return((mean(x[exp])-mean(x[cont])))}

pValues <- function (x) {return(tryCatch(t.test(x[exp],x[cont],var.equal=FALSE,paired=FALSE,conf.level = 0.95)$p.value, error=function(e) NA))}


## exp=2:3; cont=6:7

chooseID=grep("_Log",colnames(proteinRollup2))
proteinRollup3=proteinRollup2[,chooseID]

samples=colnames(proteinRollup3)
newName=gsub("EP1_meanNormalize_Log", "", samples)
newName=gsub("EP2_meanNormalize_Log", "", newName)
newName=gsub("EP3_meanNormalize_Log", "", newName)
colnames(proteinRollup3)=newName


ab=seq(1, (ncol(proteinRollup3)-2), by=3)

StatisValues <- data.frame(matrix(NA, nrow = nrow(proteinRollup3), ncol = 0));
namesStatis=c();

for(i in 1:length(ab)){
 for(j in 1:length(ab)){
     if(i!=j){	 
     name=paste(newName[ab[i]],newName[ab[j]],sep="_vs_");
     print(name);	 
	 match=grep(name,comps);
	 if(length(match)){
	   exp=ab[i]:(ab[i]+2); 
	   cont=ab[j]:(ab[j]+2);
	   ## subX=proteinRollup3[,c(exp,cont)]
       fold=apply(proteinRollup3,1,foldChange);
       pvalue=apply(proteinRollup3,1,pValues);	   
	   StatisValues=cbind(StatisValues,fold);
	   StatisValues=cbind(StatisValues,pvalue);
	   nameF=paste(name,"_fold",sep="");
	   nameP=paste(name,"_Pvalue",sep="")
	   namesStatis=c(namesStatis,nameF,nameP)
	 }    
	}
  j=j+3
  }
 i=i+3
  }


proteinRollup4=cbind(proteinRollup2,StatisValues)
name6=c(colnames(proteinRollup2),namesStatis)
colnames(proteinRollup4)=name6


write.table(proteinRollup4,"Combined_Extracellular_Master_statis_unique_pep_only_no_decoy.txt",sep="\t",col.names = NA)






