### Created by Mao Peng (m.peng@cbs.knaw.nl)on 12-01-2016
### Updated by Victoria Aguilar (v.aguilar@cbs.knaw.nl) on 18-01-2016

TTest <- function (aData,numC,numE){
cont=1:(numC);
exp=((1+numC):(numC+numE));  #### this should be careful if we have different replicate number

control=colnames(aData)[1]

sample=colnames(aData)[1+numC]
## print(c(numC,numE))


foldChange <- function (x) {return(tryCatch(mean(x[exp])/mean(x[cont]), error=function(e) NA))}

### pValues <- function (x) {return(tryCatch(t.test(x[exp],x[cont],var.equal=TRUE,paired=FALSE,conf.level = 0.95)$p.value, error=function(e) NA))}

pValues_cal<- function (x) {
 

##pValueS=tryCatch(shapiro.test(c(x[exp],x[exp]))$p.value, error=function(e) NA);
##pValueF=tryCatch(var.test(x[exp], x[cont])$p.value, error=function(e)(NA));
pValueS=tryCatch(shapiro.test(c(x[exp],x[cont]))$p.value, error=function(e) NA);
pValueF=tryCatch(var.test(x[exp],x[cont])$p.value, error=function(e) NA);
if(is.na(pValueS)){pValue=NA;return(pValue); } 
if(pValueS>=0.05){ if (pValueF >=0.05 ) {return(tryCatch(t.test(x[exp],x[cont],var.equal=TRUE,paired=FALSE,conf.level = 0.95)$p.value, error=function(e) NA)); } 
else {return(tryCatch(t.test(x[exp],x[cont],var.equal=FALSE,paired=FALSE,conf.level = 0.95)$p.value, error=function(e) NA))}
}else return(tryCatch(wilcox.test(x[cont],x[exp],correction=TRUE)$p.value, error=function(e) NA));  
}

pValues_method<-function (x) { 
pValueS=tryCatch(shapiro.test(c(x[exp],x[cont]))$p.value, error=function(e) NA);
pValueF=tryCatch(var.test(x[exp],x[cont])$p.value, error=function(e) NA);
if(is.na(pValueS)){pValue=NA;meth=NA; return(NA); } 
if(pValueS>=0.05){ 
	if (pValueF >=0.05 ) {return(tryCatch(t.test(x[exp],x[cont],var.equal=TRUE,paired=FALSE,conf.level = 0.95)$method, error=function(e) NA)); } 
	else {return(tryCatch(t.test(x[exp],x[cont],var.equal=FALSE,paired=FALSE,conf.level = 0.95)$method, error=function(e) NA))}
}else return(tryCatch(wilcox.test(x[cont],x[exp],correction=TRUE)$method, error=function(e) NA));  

}



pvalue_sha=apply(aData,1,function(x){return(tryCatch(shapiro.test(c(x[exp],x[cont]))$p.value, error=function(e) NA))});
pvalue_var=apply(aData,1,function(x){return(tryCatch(var.test(x[exp], x[cont])$p.value, error=function(e) NA))});

pvalue_t=apply(aData,1,pValues_cal);
methods =apply(aData,1,pValues_method);


fold=apply(aData,1,foldChange)
meanE=apply(aData[,exp],1,mean)
sdE=apply(aData[,exp],1,sd)
meanC=apply(aData[,cont],1,mean)
sdC=apply(aData[,cont],1,sd)


colname2=paste(control, "mean",sep="_");
colname3=paste(control,"sd",sep="_");
colname4=paste(sample,"mean",sep="_");
colname5=paste(sample,"sd",sep="_");
colname6=paste(sample,"_over_",control,sep="");
colname7="pValue_ttest";
colname8="Statistics_method";
colname9="Shapiro test P-value";
colname10="F-test P-value";

result=data.frame(meanC,sdC,meanE,sdE,fold,pvalue_t, methods, pvalue_sha, pvalue_var)
##result=data.frame(meanC,sdC,meanE,sdE,fold,pvalue_sha,pvalue_var)
colnames(result)=c(colname2,colname3,colname4,colname5,colname6,colname7,colname8,colname9,colname10)
##colnames(result)=c(colname2,colname3,colname4,colname5,colname6,colname9,colname10)
return(result)



}


