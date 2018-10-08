setwd("D:\\Projects\\proteomics_EMSL\\new_table_make2");
require(xlsx);  ## file1=read.xlsx("FE_lists_add_criteria.xlsx", sheetIndex = 1)
library(readxl) 
   
read_excel_allsheets <- function(filename, tibble = FALSE) {
    # I prefer straight data.frames but if you like tidyverse tibbles (the default with read_excel)
    # then just pass tibble = TRUE
    sheets <- readxl::excel_sheets(filename)
    x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X,col_names = TRUE))
    if(!tibble) x <- lapply(x, as.data.frame)
    names(x) <- sheets
    x
}



mysheets2 <- read_excel_allsheets("Combined_Extracellular_Master_statis_unique_pep_only_no_decoy_no_contaminant_fold2.5_renamed_v1.xlsx")
comps2=names(mysheets2)


mysheets3 <- read_excel_allsheets("DGE_summary_Sara_RNAseq_significance2_v2.xlsx")
comps3=names(mysheets3)


mysheets4 <- read_excel_allsheets("CAZyme_function_list.xlsx")
comps4=names(mysheets4)


mysheets5 <- read_excel_allsheets("Dicsqu464_1_GeneCatalog_proteins_20151220_SigP.xlsx")
comps5=names(mysheets5)


##cazy="CAZyme_table.txt";
#############for(i in 2:length(comps3)){
for(i in 1:1){
genes=mysheets4[[i]]
colns=colnames(genes);
selectID=c(1,5,2,4)
subTab=genes[selectID]
subTabC=subTab;
## write.table(subTabC,file=cazy,sep="\t",col.names=NA);
}


for(i in 1:1){
genes=mysheets5[[i]]
colns=colnames(genes);
selectID=c(1,5)
subTab=genes[selectID];
subTabS=subTab;
subTabS[,1]=paste("Dicsqu464_1_PID",subTabS[,1],sep="_");
## write.table(subTabS,file=cazy,sep="\t",col.names=NA);
}












#### for comparing 3 samples

temp1=c("average.FPKM.")
temp2=c("_log2FC","_Padj","Categorisation_")
##comp1=c("average.FPKM.p1_2w","average.FPKM.p2_2w","average.FPKM.PP_2w","p1_2w_v_PP_2w_log2FC","PP_2w_v_p1_2w_log2FC","p1_2w_v_PP_2w_Padj","p1_2w_v_PP_2w_Padj","p2_2w_v_PP_2w_log2FC","PP_2w_v_p2_2w_log2FC","p2_2w_v_PP_2w_Padj","Categorisation_p1_2w_v_PP_2w","Categorisation_p2_2w_v_PP_2w")
## cb=c("p1_2w","p2_2w","PP_2w")     ############ keys to change as the comparison you want
## cb=c("p1_2w","p2_2w","P1P2_2w")     ############ keys to change as the comparison you want
## cb=c("f1_2w","p2_2w","F1P2_2w")     ############ keys to change as the comparison you want
cb=c("PP_2w","FF_2w","F1P2_2w")     ############ keys to change as the comparison you want


cbb=paste(cb,collapse="_");
comp3=paste(cbb,"_transcriptome.txt",sep="");
cb2=c(paste(cb[1],cb[3],sep="_v_"),paste(cb[3],cb[1],sep="_v_"),paste(cb[2],cb[3],sep="_v_"),paste(cb[3],cb[2],sep="_v_"))
comp=c(paste(temp1,cb,sep=""),paste(cb2,temp2[1],sep=""),paste(cb2,temp2[2],sep=""),paste(temp2[3],cb2,sep=""))


#############for(i in 2:length(comps3)){
for(i in 1:1){
genes=mysheets3[[i]]
colns=colnames(genes);
selectID=c(1)
for(i in 1:length(comp)){
 a=grep(comp[i],colns)
 if(length(a)){
 selectID=c(selectID,a)
 }
}
## print(colns)
print(selectID)

subTabT=genes[selectID]
## write.table(subTabT,file=comp3,sep="\t",col.names=NA);
}


samples=colnames(subTabT)
newName=gsub("_log2FC", "_Fold_change", samples)
newName=gsub("_v_", "_vs_", newName)
colnames(subTabT)=newName


a=grep("_log2FC",samples)

check1=paste(cb[1],cb[3],sep="_vs_");
check2=paste(cb[2],cb[3],sep="_vs_");
replace1=paste(cb[3],cb[1],sep="_vs_");
replace2=paste(cb[3],cb[2],sep="_vs_");


for(i in 1:length(a)){
fold=subTabT[,a[i]]
fold=2^(as.numeric(fold));
for(j in 1:length(fold)){
    if(!(is.na(fold[j])) & fold[j]<1){fold[j]=-(1/fold[j])}  
  }
if(length(grep(replace1,newName[a[i]]))){fold=-fold}; 
if(length(grep(replace2,newName[a[i]]))){fold=-fold}; 
 subTabT[,a[i]]=fold
}

newName=gsub(replace1, check1, newName)
newName=gsub(replace2, check2, newName)
colnames(subTabT)=newName





temp1=c("EP1_meanNormalize\\b","EP2_meanNormalize\\b","EP3_meanNormalize\\b")
temp2=c("_fold","_Pvalue","_status")

comp3=paste(cbb,"_proteomics.txt",sep="");
cb2=c(paste(cb[1],cb[3],sep="_vs_"),paste(cb[3],cb[1],sep="_vs_"),paste(cb[2],cb[3],sep="_vs_"),paste(cb[3],cb[2],sep="_vs_"))
comp=c(paste(cb[1],temp1,sep=""),paste(cb[2],temp1,sep=""),paste(cb[3],temp1,sep=""),paste(cb2,temp2[1],sep=""),paste(cb2,temp2[2],sep=""),paste(cb2,temp2[3],sep=""))

#############for(i in 2:length(comps2)){
for(i in 1:1){
genes=mysheets2[[i]]
colns=colnames(genes);
selectID=c(2)
for(i in 1:length(comp)){
 a=grep(comp[i],colns)
 if(length(a)){
 selectID=c(selectID,a)
 }
}
## print(colns)
print(selectID)

subTab=genes[selectID]
average1=apply(subTab[,2:4],1,mean);
average2=apply(subTab[,5:7],1,mean);
average3=apply(subTab[,8:10],1,mean);
newNames=gsub("2wEP1_meanNormalize","2w_average_abundance",colnames(subTab)[c(2,5,8)],)
name2=colnames(subTab)

subTab$average1=average1;
subTab$average2=average2;
subTab$average3=average3;
colnames(subTab)=c(name2,newNames)
subTabP=subTab;
## write.table(subTabP,file=comp3,sep="\t",col.names=NA);
}




samples=colnames(subTabP)
newName=gsub("_fold", "_Fold_Change", samples)
newName=gsub("_v_", "_vs_", newName)
newName=gsub("_status", "_categorisation", newName)
colnames(subTabP)=newName

a=grep("_Fold_Change",newName)

check1=paste(cb[1],cb[3],sep="_vs_");
check2=paste(cb[2],cb[3],sep="_vs_");
replace1=paste(cb[3],cb[1],sep="_vs_");
replace2=paste(cb[3],cb[2],sep="_vs_");


for(i in 1:length(a)){
fold=subTabP[,a[i]]
fold=as.numeric(fold);
for(j in 1:length(fold)){
    if(!(is.na(fold[j])) & fold[j]<1){fold[j]=-(1/fold[j])}  
  }
if(length(grep(replace1,newName[a[i]]))){fold=-fold}; 
if(length(grep(replace2,newName[a[i]]))){fold=-fold}; 
 subTabP[,a[i]]=fold
}

newName=gsub(replace1, check1, newName)
newName=gsub(replace2, check2, newName)
colnames(subTabP)=newName

subTabP2=subTabP[,c(1,17,18,19,11,12,13,14,15,16)];
comp4=paste(cbb,"_comparison",sep="");
comb=merge(subTabT,subTabP2,by.x="Dicsqu464_1_PID_format",by.y="Dicsqu464_1_PID_format",all=TRUE);
comb=merge(subTabC,comb,by.x="Dicsqu464_1 proteins",by.y="Dicsqu464_1_PID_format");
comb=merge(subTabS,comb,by.x="proteinid",by.y="Dicsqu464_1 proteins",all.y=TRUE);


selectID=c();
for(i in 1:nrow(comb)){
 a=grep("Greater_in_",comb[i,])
 if(length(a)){
 selectID=c(selectID,i)
}
}
comb2=comb[selectID,]
## write.table(comb2,file=comp4,sep="\t",col.names=NA);	


write.xlsx(x = comb2, file = "combined_transcriptomics_proteomics_data.xlsx",
        sheetName =comp4, row.names = FALSE, append=TRUE)
		
		
		


		
		

		
		
		
		




		

		
		
#### for comapring 2 samples

temp1=c("average.FPKM.")
temp2=c("_log2FC","_Padj","Categorisation_")
##comp1=c("average.FPKM.p1_2w","average.FPKM.p2_2w","average.FPKM.PP_2w","p1_2w_v_PP_2w_log2FC","PP_2w_v_p1_2w_log2FC","p1_2w_v_PP_2w_Padj","p1_2w_v_PP_2w_Padj","p2_2w_v_PP_2w_log2FC","PP_2w_v_p2_2w_log2FC","p2_2w_v_PP_2w_Padj","Categorisation_p1_2w_v_PP_2w","Categorisation_p2_2w_v_PP_2w")
## cb=c("FF_2w","f1_2w")     ############ keys to change as the comparison you want
## cb=c("PP_2w","FF_2w")     ############ keys to change as the comparison you want
cb=c("PP_2w","P1P2_2w")     ############ keys to change as the comparison you want


cbb=paste(cb,collapse="_");
comp3=paste(cbb,"_transcriptome.txt",sep="");
cb2=c(paste(cb[1],cb[2],sep="_v_"),paste(cb[2],cb[1],sep="_v_"))
comp=c(paste(temp1,cb,sep=""),paste(cb2,temp2[1],sep=""),paste(cb2,temp2[2],sep=""),paste(temp2[3],cb2,sep=""))



#############for(i in 2:length(comps3)){
for(i in 1:1){
genes=mysheets3[[i]]
colns=colnames(genes);
selectID=c(1)
for(i in 1:length(comp)){
 a=grep(comp[i],colns)
 if(length(a)){
 selectID=c(selectID,a)
 }
}
## print(colns)
print(selectID)

subTabT=genes[selectID]
## write.table(subTabT,file=comp3,sep="\t",col.names=NA);
}


samples=colnames(subTabT)
newName=gsub("_log2FC", "_Fold_change", samples)
newName=gsub("_v_", "_vs_", newName)
colnames(subTabT)=newName


a=grep("_log2FC",samples)

check1=paste(cb[1],cb[2],sep="_vs_");
replace1=paste(cb[2],cb[1],sep="_vs_");


for(i in 1:length(a)){
fold=subTabT[,a[i]]
fold=2^(as.numeric(fold));
for(j in 1:length(fold)){
    if(!(is.na(fold[j])) & fold[j]<1){fold[j]=-(1/fold[j])}  
  }
if(length(grep(replace1,newName[a[i]]))){fold=-fold};  
 subTabT[,a[i]]=fold
}

newName=gsub(replace1, check1, newName)
colnames(subTabT)=newName




temp1=c("EP1_meanNormalize\\b","EP2_meanNormalize\\b","EP3_meanNormalize\\b")
temp2=c("_fold","_Pvalue","_status")
##comp1=c("average.FPKM.p1_2w","average.FPKM.p2_2w","average.FPKM.PP_2w","p1_2w_v_PP_2w_log2FC","PP_2w_v_p1_2w_log2FC","p1_2w_v_PP_2w_Padj","p1_2w_v_PP_2w_Padj","p2_2w_v_PP_2w_log2FC","PP_2w_v_p2_2w_log2FC","p2_2w_v_PP_2w_Padj","Categorisation_p1_2w_v_PP_2w","Categorisation_p2_2w_v_PP_2w")

comp3=paste(cbb,"_proteomics.txt",sep="");
cb2=c(paste(cb[1],cb[2],sep="_vs_"),paste(cb[2],cb[1],sep="_vs_"))
comp=c(paste(cb[1],temp1,sep=""),paste(cb[2],temp1,sep=""),paste(cb2,temp2[1],sep=""),paste(cb2,temp2[2],sep=""),paste(cb2,temp2[3],sep=""))

#############for(i in 2:length(comps2)){
for(i in 1:1){
genes=mysheets2[[i]]
colns=colnames(genes);
selectID=c(2)
for(i in 1:length(comp)){
 a=grep(comp[i],colns)
 if(length(a)){
 selectID=c(selectID,a)
 }
}
## print(colns)
print(selectID)

subTab=genes[selectID]
average1=apply(subTab[,2:4],1,mean);
average2=apply(subTab[,5:7],1,mean);
newNames=gsub("2wEP1_meanNormalize","2w_average_abundance",colnames(subTab)[c(2,5)],)
name2=colnames(subTab)

subTab$average1=average1;
subTab$average2=average2;

colnames(subTab)=c(name2,newNames)
subTabP=subTab;
## write.table(subTabP,file=comp3,sep="\t",col.names=NA);
}



samples=colnames(subTabP)
newName=gsub("_fold", "_Fold_Change", samples)

newName=gsub("_status", "_categorisation", newName)

colnames(subTabP)=newName

a=grep("_Fold_Change",newName)

check1=paste(cb[1],cb[2],sep="_vs_");
replace1=paste(cb[2],cb[1],sep="_vs_");


for(i in 1:length(a)){
fold=subTabP[,a[i]]
fold=as.numeric(fold);
for(j in 1:length(fold)){
    if(!(is.na(fold[j])) & fold[j]<1){fold[j]=-(1/fold[j])}  
  }
if(length(grep(replace1,newName[a[i]]))){fold=-fold};  
 subTabP[,a[i]]=fold
}

newName=gsub(replace1, check1, newName)
colnames(subTabP)=newName


subTabP2=subTabP[,c(1,11,12,8,9,10)]
comp4=paste(cbb,"_comparison",sep="");
comb=merge(subTabT,subTabP2,by.x="Dicsqu464_1_PID_format",by.y="Dicsqu464_1_PID_format",all=TRUE);
comb=merge(subTabC,comb,by.x="Dicsqu464_1 proteins",by.y="Dicsqu464_1_PID_format");
comb=merge(subTabS,comb,by.x="proteinid",by.y="Dicsqu464_1 proteins",,all.y=TRUE);


selectID=c();
for(i in 1:nrow(comb)){
 a=grep("Greater_in_",comb[i,])
 if(length(a)){
 selectID=c(selectID,i)
}
}
comb2=comb[selectID,]

write.xlsx(x = comb2, file = "combined_transcriptomics_proteomics_data.xlsx",
        sheetName =comp4, row.names = FALSE, append=TRUE)


