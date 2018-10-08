setwd("D:\\Projects\\proteomics_EMSL\\Master_sheets\\Extracellular");

data=read.table("Combined_Extracellular_Master_Sheet_ZeroFill.txt",sep="\t",head=T);

colnames(data)
x=as.matrix(data[,5:ncol(data)])   ### be sure which colum start to have the real instensity values, here is colum-5  
adjust=mean(colMeans(x));

x2=x;
for(i in 1:ncol(x)){
x2[,i]=x[,i]*(adjust/mean(x[,i]))
}

x3=log2(x2);
x3[is.infinite(x3)] <- NA

data[,5:ncol(data)]=x3

write.table(data,"Combined_Extracellular_Master_Sheet_pepNormalized.txt",sep="\t",col.names = NA)


x3[is.na(x3)] <- 0
samples=c();
averageM=matrix(nrow=(dim(x3))[1],ncol=ncol(x3)/3);
## meanPep <- function (x) {a=is.na(x); if(sum(a)>=2){return(NA)}; if(sum(a)<2){return(mean(x[!a]))} }
meanPep <- function (x) {a=(x==0); if(sum(a)>=2){return(0)}; if(sum(a)<2){return(mean(x[!a]))} }

j=0;
for(k in 1:ncol(x3)){
if(k%%3==1){
  j=j+1;
  samples=c(samples,(colnames(x3))[k])
  b=x3[,k:(k+2)];
  meanExp=apply(b,1,meanPep);
  averageM[,j]=meanExp;
}
}



samples=gsub("EP1", "EP_average", samples)
colnames(averageM)=samples

deLogAverage=2^averageM;
averageM2=cbind(data[,1:4],averageM)
write.table(averageM2,"Combined_Extracellular_Master_Sheet_averagePep.txt",sep="\t",col.names = NA)

deLogAverage=2^averageM;
deLogAverage[deLogAverage==1] <- 0
deLogAverage2=cbind(data[,1:4],deLogAverage)
write.table(deLogAverage2,"Combined_Extracellular_Master_Sheet_DelogAveragePep.txt",sep="\t",col.names = NA)


deLogAverage3=deLogAverage2[,c(3,5:ncol(deLogAverage2))]
proteinRollup= aggregate(. ~ Reference, data=deLogAverage3, FUN=sum)




comparison=read.table("comparisons.txt",sep="\t",head=T);
comps=comparison$comparisons


## samples=colnames(proteinRollup)
## EP2=grep("2EP",samples); EP4=grep("4EP",samples);

foldChange <- function (x) { if(x[1]>0 & x[2]>0){return(x[1]-x[2])}else{return(NA)} }


proteinRollup2=proteinRollup;
proteinRollup2[proteinRollup2==0] <- 1
proteinRollup2[,2:ncol(proteinRollup2)]=log2(proteinRollup2[,2:ncol(proteinRollup2)]);

samples=colnames(proteinRollup2)
newName=gsub("EP_average", "", samples)
colnames(proteinRollup2)=newName


for(i in 2:(ncol(proteinRollup2)-1)){
 for(j in 3:(ncol(proteinRollup2))){
     if(i!=j){
     name=paste(newName[i],newName[j],sep="_vs_")	 
	 match=grep(name,comps);
	 if(length(match)){
	   subX=proteinRollup2[,c(i,j)]
       fold=apply(subX,1,foldChange);	     
	   proteinRollup2=cbind(proteinRollup2,fold);
	   newName=c(newName,name)
	 }
    }
  }
}
colnames(proteinRollup2)=newName



write.table(proteinRollup2,"Combined_Extracellular_Master_foldChange.txt",sep="\t",col.names = NA)






proteinRollup3=proteinRollup2;
compID=grep("_vs_",colnames(proteinRollup3));

####
Names=colnames(proteinRollup3);
StatisName=c();
StatisValues <- data.frame(matrix(NA, nrow = nrow(proteinRollup3), ncol = 0));


#### plot the distribution within a pdf file
pdf("Fold_change_Normalization_check_Extracellular.pdf")
## dev.off()

lod1=0;
lod2=0;
sigma2=0;
Sample1="s1";
Sample2="s2";
Significance <- function (x) { if(((x[1])==0) & ((x[2])==0)){return("Not present")}
                              else if(((x[2])==0) & ((x[1]>0))){   if(x[1]>lod1){return(paste("Only_in",Sample1,sep="_")) }else{return("Not Significant")}      }
                              else if(((x[1])==0) & ((x[2]>0))){   if(x[2]>lod2){return(paste("Only_in",Sample2,sep="_"))}else{return("Not Significant")}  }
							  else if(((x[3]<0)) & (abs(x[3])>sigma2)){   {return(paste("Greater_in",Sample2," Fold_change",round(2^(abs(x[3])),digits=2),sep="_"))}  }
							  else if(((x[3]>0)) & (x[3]>sigma2)){   {return(paste("Greater_in",Sample1," Fold_change",round(2^((x[3])),digits=2), sep="_"))}  }
                              else {return("Not Significant") }
							   }
##### x[1] x[2] - real intensity;||  lod1, lod2 -  lower detected intensity (Mean-2sigma);  ||   x[3]  - Normalized fold change  || sigma2-  two sigma of fold change




for(i in 1:length(compID)){
k=compID[i];
compName=colnames(proteinRollup2)[k];
tmp=strsplit(compName, "_vs_")[[1]]
Sample1=tmp[1];
Sample2=tmp[2];
id1=grep(Sample1,colnames(proteinRollup2))[1];
id2=grep(Sample2,colnames(proteinRollup2))[1];

temp=proteinRollup2[,id1];
temp=temp[(temp)>0];
lod1=mean(temp)-2*(sd(temp))

temp=proteinRollup2[,id2];
## temp=temp[!(is.na(temp))];
temp=temp[(temp)>0];
lod2=mean(temp)-2*(sd(temp))

### lod1=100; lod2=100;   ### test the globlal vairable works or not
print(c(Sample1,Sample2));   print(c(lod1,lod2)); 

normalizedName=paste("Normalized",compName,sep="_");
Names=c(Names,normalizedName);
statis=paste(compName,"status",sep=" ");
StatisName=c(StatisName,statis);


test=proteinRollup2[,k];
test2=test[!(is.na(test))];
myhist=hist(test2, breaks=seq(min(test2)-0.5,max(test2)+0.5,0.5),xlim=c(-10,10),plot=FALSE)

multiplier <- myhist$counts / myhist$density
mydensity <- density(test2)
mydensity$y <- mydensity$y * multiplier[1]

adjust=median(test2)
test3=test2-adjust
sigma2=2*(sd(test3))

test4=test-adjust
proteinRollup3=cbind(proteinRollup3,test4)
myhist2=hist(test3, breaks=seq(min(test3)-0.5,max(test3)+0.5,0.5),xlim=c(-10,10),plot=FALSE)

multiplier2 <- myhist2$counts / myhist2$density
mydensity2 <- density(test3)
mydensity2$y <- mydensity2$y * multiplier2[1]

plot(myhist,xlim=c(-10,10), main=compName)
lines(mydensity,col="red",lwd=2)

lines(mydensity2,col="green",lwd=2)
### abline(v =0,col="orange",lwd=2)
abline(v =median(test2),col="red",lty=2,lwd=1.5)
abline(v =median(test3),col="green",lty=2,lwd=1.5)

adjust2=-round(adjust,digits=3)
adjust3=paste("Adjust_value",adjust2,sep=": ")
legend("topleft",c("Before Normalize","After Normalize",adjust3),lwd=2,lty=c(1,1,3),col=c("red","green","blue"),bty="n")


subY=proteinRollup2[,c(id1,id2)]
subY=cbind(subY,test4)

sig=apply(subY,1,Significance);
StatisValues=cbind(StatisValues,sig)
}




dev.off()

Names=c(Names,StatisName);
proteinRollup4=cbind(proteinRollup3,StatisValues)
colnames(proteinRollup4)=Names
write.table(proteinRollup4,"Combined_Extracellular_Master_foldChange_statis.txt",sep="\t",col.names = NA)


