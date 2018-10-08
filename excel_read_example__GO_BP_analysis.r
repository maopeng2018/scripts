setwd("D:\\test\\GO_analysis");

###  GO annotation for the whole genome
GOs=read.table("Aspni_NRRL3_GO_CAZyme_KEGG.txt",sep="\t",head=T,stringsAsFactors = FALSE,quote="");
DsGO=GOs$GO_BP
DsGO=strsplit(DsGO, ";")
DsGO=unlist(DsGO)
DsGO=gsub("\"","",DsGO)
GOs$X=gsub("jgi.p\\|Aspni_NRRL3_1\\|","",GOs$X)


## require(xlsx);  ## file1=read.xlsx("FE_lists_add_criteria.xlsx", sheetIndex = 1)
library(readxl)    
read_excel_allsheets <- function(filename, tibble = FALSE) {
    # I prefer straight data.frames but if you like tidyverse tibbles (the default with read_excel)
    # then just pass tibble = TRUE
    sheets <- readxl::excel_sheets(filename)
    x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X,col_names = FALSE))
    if(!tibble) x <- lapply(x, as.data.frame)
    names(x) <- sheets
    x
}

## mysheets <- read_excel_allsheets("FE_lists_add_criteria.xlsx")
mysheets <- read_excel_allsheets("FE_lists_DevIN_2FC_P0.05.xlsx")

comps=names(mysheets)
for(i in 1:length(comps)){
    genes=mysheets[[i]]
	genes[,1]=gsub("NRRL3_0000","",genes[,1])   #### keep the replace order from _000000  to _0 , otherwise there will be problem 
	genes[,1]=gsub("NRRL3_000","",genes[,1])
	genes[,1]=gsub("NRRL3_00","",genes[,1])
	genes[,1]=gsub("NRRL3_0","",genes[,1])
    genes[,1]=gsub("NRRL3_","",genes[,1])
    print(comps[i]); print( dim(genes))

goAdd=merge(genes,GOs,by.x="X__1",by.y="X")	
 print(goAdd[1:2,])
GOup=goAdd$GO_BP
GOup=strsplit(GOup, ";")
GOup=unlist(GOup)
GOup=gsub("\"","",GOup) 


upInfor=goAdd[,c("X__1","GO_BP")]
 

comp3=paste(comps[i],"_GO_BP_enrichments.txt",sep="_")


GOup_count=sort(table(GOup))
DsGO_count=sort(table(DsGO))

table1=cbind(terms=names(DsGO_count), count1=unname(unlist(DsGO_count)))
table2=cbind(terms=names(GOup_count), count2=unname(unlist(GOup_count)))


##### for up-regulated gene enrichment analysis
compare=merge(table1,table2,by="terms")
compare2=compare[,2:3]
compare2=as.matrix(compare2)
compare2=apply(compare2, 2, as.numeric)
row.names(compare2)=compare[,1]
### probToDrawWhite=phyper(drawnWhite,TotalWhiteBall,TotalBlackBall,drawnBall)

DrawnBall=dim(genes)[1]
wtDrawnID=2
wtTotalID=1
TotalBall=(dim(GOs))[1]


pValues_cal<- function (x) { 
blackBall=TotalBall-as.numeric(x[wtTotalID]);
##pval<-return(tryCatch(phyper(x[wtDrawnID],x[wtTotalID],(TotalBall-x[wtTotalID]),DrawnBall),error=function(e)NA))
return(1-phyper(x[wtDrawnID],x[wtTotalID],blackBall,DrawnBall))
}

foldChange <- function (x) {return((x[wtDrawnID]/DrawnBall)/(x[wtTotalID]/TotalBall))}


pvalues=apply(compare2,1,pValues_cal);
folds=apply(compare2,1,foldChange);

result=data.frame(compare2,folds,pvalues)
colnames(result)=c("genes_in_wholeGenome","genes_in_significantChange_list"，"enriched_fold","p-value")
result2=result[result[,2]>=1,]

genes=c();
names=c();
nameF=paste(comps[i],"_genes",sep="");

for(d in 1:nrow(result2)){
gene="/";
go=rownames(result2)
ids=grep(go[d],upInfor[,2],fixed=TRUE);
if(length(ids)){
gene=paste(upInfor[ids,1],collapse=";")
}
genes=c(genes,gene)
}
StatisValues=cbind(result2,genes);
namesStatis=c(colnames(result2),nameF)
colnames(StatisValues)=namesStatis

write.table(StatisValues,file=comp3,sep="\t",col.names=NA);

 }

