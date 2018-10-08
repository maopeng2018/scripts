### Created by Mao Peng (m.peng@cbs.knaw.nl)on 12-01-2016



##pValues <- function (x) {return(tryCatch(t.test(x[exp],x[cont],var.equal=TRUE,paired=FALSE,conf.level=0.95)$p.value, error=function(e) NULL))}






TTest <- function (aData,numC,numE){
cont=1:(numC);
exp=((1+numC):(numC+numE));  #### this should be careful if we have different replicate number

control=colnames(aData)[1]

sample=colnames(aData)[1+numC]
## print(c(numC,numE))


foldChange <- function (x) {return(tryCatch(mean(x[exp])/mean(x[cont]), error=function(e) NA))}

### pValues <- function (x) {return(tryCatch(t.test(x[exp],x[cont],var.equal=TRUE,paired=FALSE,conf.level = 0.95)$p.value, error=function(e) NA))}

pValues_cal<- function (x) { 
pval<-return(tryCatch(t.test(x[exp],x[cont],var.equal=TRUE,paired=FALSE,conf.level = 0.95)$p.value,error=function(e)NA))

}

pvalues=apply(aData,1,pValues_cal);

fold=apply(aData,1,foldChange)
meanE=apply(aData[,exp],1,mean)
sdE=apply(aData[,exp],1,sd)
 meanC=apply(aData[,cont],1,mean)
 sdC=apply(aData[,cont],1,sd)


colname2=paste(control,"mean",sep="_");
colname3=paste(control,"sd",sep="_");
colname4=paste(sample,"mean",sep="_")
colname5=paste(sample,"sd",sep="_")
colname6=paste(sample,"_over_",control,sep="");
colname7="pValue_ttest";

result=data.frame(meanC,sdC,meanE,sdE,fold,pvalues)
colnames(result)=c(colname2,colname3,colname4,colname5,colname6,colname7)
return(result)



}


