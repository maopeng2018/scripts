setwd("D:\\Projects\\Kristiina_D_squalens\\Mao_analysis\\co_culture\\HTseq_result\\DEseq2")
data=read.table("All_counts_with_annotation_statis.txt",sep="\t",head=T,stringsAsFactors = FALSE,quote="");

dim(data)
ids=data[,1]
Ds=data[grep("Dicsqu",ids),]
Fp=data[grep("Fompi",ids),]

comp="Ds2w_VS_Ds3w"
comp2=paste(comp,"Pvalue",sep="_")
comp3=paste(comp,"_upRegulatedGene_GO_enrichments.txt",sep="_")
comp4=paste(comp,"_downRegulatedGene_GO_enrichments.txt",sep="_")
choseColum=grep(comp2,colnames(Ds))

values=Ds[,choseColum]
NAs=which(is.na(values), TRUE)
zero=grep("\\\\",values)
Ds2=Ds[-(c(NAs,zero)),]

values=Ds2[,choseColum]
values=as.numeric(as.character(values))
select=Ds2[values<0.01,]
### write.table(select,file="p0001.txt",sep="\t",col.names=NA);

choseColum2=choseColum-1
up_Total=sum(select[,choseColum2]>1)
down_Total=sum(select[,choseColum2]<1)

cazyID=grep("CAZyme",colnames(select))
zero=grep("\\\\",select[,cazyID])
GT=grep("GlycosylTransferase",select[,cazyID])

cazymes=select[-c(zero,GT),]
cazymes[,cazyID]=gsub("GlycosylTransferase Family ","GT",cazymes[,cazyID])
cazymes[,cazyID]=gsub("Glycosyltransferase Family ","GT",cazymes[,cazyID])
cazymes[,cazyID]=gsub("Carbohydrate-Binding Module Family ","CBM",cazymes[,cazyID])
cazymes[,cazyID]=gsub("Auxiliary Activity Family ","AA",cazymes[,cazyID])
cazymes[,cazyID]=gsub("Glycoside Hydrolase Family ","GH",cazymes[,cazyID])
cazymes[,cazyID]=gsub(" / Subf ","_",cazymes[,cazyID])
cazymes[,cazyID]=gsub("Polysaccharide Lyase Family ","PL",cazymes[,cazyID])
cazymes[,cazyID]=gsub("Carbohydrate Esterase Family ","CE",cazymes[,cazyID])
cazymes[,cazyID]=gsub("Distantly related to plant expansins","expansins",cazymes[,cazyID])


cazymeUp=cazymes[cazymes[,choseColum2]>1,cazyID]
cazymeDown=cazymes[cazymes[,choseColum2]<1,cazyID]

intersect(cazymeUp,cazymeDown); 
setdiff(cazymeUp,cazymeDown);
setdiff(cazymeDown,cazymeUp);
sort(table(cazymeUp))


library("gplots")
##library("heatmap.plus")

values2=log2(cazymes[,7:12])
values2=as.matrix(values2)
names=paste(cazymes$SequenceID,cazymes$CAZyme,sep = "_");
names=gsub("Dicsqu464_1_","Ds_",names)

row.names(values2)=names
	  
heatmap.2(values2,
 margin=c(15,15),
		     col=bluered(75),
             Colv=NA,
			 key=F,
			 trace="none",
			 cexCol=1.2,
			 cexRow=0.5,
             ylab="Differential expressed Cazymes",
) 




dim(select)
selectUp=select[select[,choseColum2]>1,]
selectDown=select[select[,choseColum2]<1,]

GOup=selectUp$GO_MF
GOdown=selectDown$GO_MF

GOup=strsplit(GOup, ";")
GOup=unlist(GOup)
GOup=gsub("\"","",GOup)

GOdown=strsplit(GOdown, ";")
GOdown=unlist(GOdown)
GOdown=gsub("\"","",GOdown)


DsGO=Ds$GO_MF
DsGO=strsplit(DsGO, ";")
DsGO=unlist(DsGO)
DsGO=gsub("\"","",DsGO)

GOup_count=sort(table(GOup))
GOdown_count=sort(table(GOdown))
DsGO_count=sort(table(DsGO))

table1=cbind(terms=names(DsGO_count), count1=unname(unlist(DsGO_count)))
table2=cbind(terms=names(GOup_count), count2=unname(unlist(GOup_count)))
table3=cbind(terms=names(GOdown_count), count3=unname(unlist(GOdown_count)))

##### for up-regulated gene enrichment analysis

compare=merge(table1,table2,by="terms")
compare2=compare[,2:3]
compare2=as.matrix(compare2)
compare2=apply(compare2, 2, as.numeric)
row.names(compare2)=compare[,1]
### probToDrawWhite=phyper(drawnWhite,TotalWhiteBall,TotalBlackBall,drawnBall)

DrawnBall=(dim(select))[1]
wtDrawnID=2
wtTotalID=1
TotalBall=(dim(Ds))[1]

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
result2=result[result[,2]>=5,]
write.table(result2,file=comp3,sep="\t",col.names=NA);


##### for down-regulated gene enrichment analysis


compare=merge(table1,table3,by="terms")
compare2=compare[,2:3]
compare2=as.matrix(compare2)
compare2=apply(compare2, 2, as.numeric)
row.names(compare2)=compare[,1]
### probToDrawWhite=phyper(drawnWhite,TotalWhiteBall,TotalBlackBall,drawnBall)

DrawnBall=(dim(select))[1]
wtDrawnID=2
wtTotalID=1
TotalBall=(dim(Ds))[1]

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
result2=result[result[,2]>=5,]
write.table(result2,file=comp4,sep="\t",col.names=NA);


