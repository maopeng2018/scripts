setwd("D:\\Projects\\proteomics_EMSL\\Master_sheets\\Extracellular_ttest_uniquePep_without_contaminant2\\test_RNA-seq");

data=read.table("DGE_summary_Sara_RNAseq.txt",sep="\t",head=T);



Sample1="s1";
Sample2="s2";

compID=grep("_v_",colnames(data));


StatisValues <- data.frame(matrix(NA, nrow = nrow(data), ncol = 0));
namesStatis=c();

ab=seq(1, (length(compID)-1), by=3)
FC=log2(2.5)  ### fold change cutoff


Significance <- function (x) {
                              mode(x) <- "numeric"
                              if( (is.na(x[1])) |  (is.na(x[2])) ){ return("Lowly expressed in both conditions") } 
                              else if( (x[1]>=FC) & ( x[2]<0.05) & (x[3]>=10) ){ return(paste("Greater_in",Sample1,sep="_")) }
                              else if( (x[1]<=(-FC)) & (x[2]<0.05) & (x[4]>=10) ){   return(paste("Greater_in",Sample2,sep="_"))       }
                              else if( (x[2]>=0.05) & ((x[3]>=10) | (x[4]>=10))  ){   return( paste("No significant difference_") )    }
                              else if( ( x[1]>(-FC)) & (x[1]<FC) & ((x[3]>=10) | (x[4]>=10))  ){   return( paste("No significant difference") )    }                              
                              else if( (x[3]<10) & (x[4]<10)  ){   return( paste("Lowly expressed in both conditions") )     }                              
							  else{	
                                      print(c(x[1],x[2],x[3],x[4],(-FC)));							  
							          return("unknown_situation");                                   
                                   }
							   }

for(i in 1:length(ab)){

       foldID=compID[ab[i]];
       pvalueID=compID[ab[i]+1];
	   name=colnames(data)[foldID]
	   name=gsub("_log2FC", "", name)
	   tmp=strsplit(name, "_v_")[[1]]
       Sample1=tmp[1];
       Sample2=tmp[2];
       id1=grep(Sample1,colnames(data))[1];
       id2=grep(Sample2,colnames(data))[1];
	    ##print(c(name,Sample1,Sample2,colnames(data)[foldID],colnames(data)[pvalueID]));	   
	   subY=data[,c(foldID,pvalueID,id1,id2)]
	   sig=apply(subY,1,Significance);
	   
	   StatisValues=cbind(StatisValues,sig);
      
	   nameSig=paste("Categorisation",name,sep="_")	   
	   namesStatis=c(namesStatis,nameSig)

  }


data2=cbind(data,StatisValues)
name6=c(colnames(data),namesStatis)
colnames(data2)=name6  


write.table(data2,"DGE_summary_Sara_RNAseq_significance2.txt",sep="\t",col.names = NA)






