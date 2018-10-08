setwd("D:\\Projects\\Evy_mycelicphlicia\\Daniel_data\\statis_analysis\\target_analysis")
getwd();


data2=read.table("comparison_all.txt",sep="\t",head=T)
cazyList=read.table("CAZY_Spoth2_Myceliophthora_thermophila.txt",sep="\t",head=T)
transporter=read.table("Transporters List JGI.txt",sep="\t",head=T)
metabolic=read.table("metabolic_Targeted_analysis_gene_lists.txt",sep="\t",head=T)


colnames(data2)
library("pheatmap")
library("gplots")
library("heatmap.plus")

ids=c("wtg2_over_mtg2","pVal.3","mtxyl2_over_wtxyl2","pVal.2","mtara2_over_wtara2","pVal.26","wtax2_over_mtax2","pVal.41","mtax8_over_wtax8","pVal.36");

mtVSwt=data2[,ids]
rownames(mtVSwt)=data2$gene.id
print((rownames(mtVSwt))[1:3])




exp_cazy=merge(mtVSwt,cazyList,by.x="row.names",by.y="Protein_Id")


RPKM=read.table("expression_RPKM.txt",sep="\t",head=T)
cazys=merge(exp_cazy,RPKM,by.x="Row.names",by.y="GeneID")


write.table(cazys,file="All_cazy_expression_profiles.txt",sep="\t",col.names=NA)




