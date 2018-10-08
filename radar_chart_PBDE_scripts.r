setwd("D:\\Projects\\creA_project\\analysis\\visualize\\CAZyme")
FPKMstatis <- read.table("RPKMs_significant_genes_fold2p001.txt",sep="\t",head=T, na.strings="NA", dec=".", strip.white=TRUE)
cazyme=read.table("A_niger_CAZymes_substrate_Paul.txt",head=T,sep="\t", na.strings="NA", dec=".", strip.white=TRUE)


combine <- merge(FPKMstatis,cazyme,by.x="JGI_id",by.y="Gene_ID") ;


data=combine

Preculture=apply(data[,3:4],1,mean);
mPreculture=apply(data[,5:6],1,mean);
W_4h=apply(data[,7:8],1,mean);
mW_4h=apply(data[,9:11],1,mean);
W_24h=apply(data[,12:14],1,mean);
mW_24h=apply(data[,15:17],1,mean);
W_48h=apply(data[,18:20],1,mean);
mW_48h=apply(data[,21:22],1,mean);
SB_4h=apply(data[,23:25],1,mean);
mSB_4h=apply(data[,26:28],1,mean);
SB_24h=apply(data[,29:31],1,mean);
mSB_24h=apply(data[,32:34],1,mean);
SB_48h=apply(data[,35:36],1,mean);
mSB_48h=apply(data[,37:39],1,mean);




## m=cbind(Preculture,mPreculture,C_4h,C_24h,C_48h,mC_4h,mC_24h,mC_48h,S_4h,S_24h,S_48h,mS_4h,mS_24h,mS_48h,SB_4h,SB_24h,SB_48h,mSB_4h,mSB_24h,mSB_48h,W_4h,W_24h,W_48h,mW_4h,mW_24h,mW_48h)
m=cbind(Preculture,mPreculture,W_4h,mW_4h,W_24h,mW_24h,W_48h,mW_48h,SB_4h,mSB_4h,SB_24h,mSB_24h,SB_48h,mSB_48h)

rownames(m)=combine$JGI_id
m_max <- (apply(m, 1, max))
m_max=unname(unlist(m_max))
## m2=m[m_max>20,]
m2=m[m_max>30,]

FPKMmax30IDs=rownames(m2)

row.names(combine)=combine$JGI_id
combine2=combine[FPKMmax30IDs,]

select=(grep("\\bno\\b",combine2$Substrate_Paul))
PBDE=combine2[-select,]

write.table(PBDE,file="RPKMs_significant_PBDE_maxFPKM30.txt",sep="\t", col.names = NA)



subPBDE=PBDE[,c(78:84,92)]
print(colnames(subPBDE))

mPreculture=sort(table(subPBDE[grep("\\bGreater_in_m",subPBDE[,1]),8]))
mW4h=sort(table(subPBDE[grep("\\bGreater_in_m",subPBDE[,2]),8]))
mW24h=sort(table(subPBDE[grep("\\bGreater_in_m",subPBDE[,3]),8]))
mW48h=sort(table(subPBDE[grep("\\bGreater_in_m",subPBDE[,4]),8]))
mSB4h=sort(table(subPBDE[grep("\\bGreater_in_m",subPBDE[,5]),8]))
mSB24h=sort(table(subPBDE[grep("\\bGreater_in_m",subPBDE[,6]),8]))
mSB48h=sort(table(subPBDE[grep("\\bGreater_in_m",subPBDE[,7]),8]))


table1=cbind(terms=names(mPreculture), count1=unname(unlist(mPreculture)))
table2=cbind(terms=names(mW4h), count1=unname(unlist(mW4h)))
table3=cbind(terms=names(mW24h), count1=unname(unlist(mW24h)))
table4=cbind(terms=names(mW48h), count1=unname(unlist(mW48h)))

table5=cbind(terms=names(mSB4h), count1=unname(unlist(mSB4h)))
table6=cbind(terms=names(mSB24h), count1=unname(unlist(mSB24h)))
table7=cbind(terms=names(mSB48h), count1=unname(unlist(mSB48h)))






compare=merge(table2,table3,by="terms")
compare=merge(compare,table4,by="terms")
compare2=compare[,2:4]
compare2=as.matrix(compare2)
compare2=apply(compare2, 2, as.numeric)
row.names(compare2)=compare[,1]
compareW=compare2[rowSums(compare2)>0,]



compare=merge(table5,table6,by="terms")
compare=merge(compare,table7,by="terms")
compare2=compare[,2:4]
compare2=as.matrix(compare2)
compare2=apply(compare2, 2, as.numeric)
row.names(compare2)=compare[,1]
compareSB=compare2[rowSums(compare2)>0,]


#### Radar plot example

library(fmsb)


compareP=t(compareW)
data=as.data.frame(compareP)

fb=data.frame(max=rep(max(data),ncol(data)),min=rep(min(data),ncol(data)))
fb=t(fb)
colnames(fb)=colnames(data)
data=rbind(fb,data)

colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.7,0.5,0.1,0.9) )
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.7,0.5,0.1,0.4) )
radarchart( data  , axistype=1 , 
    #custom polygon
    pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
    #custom the grid
    cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,20,5), cglwd=0.8,
    #custom labels
    vlcex=0.8 
    )
	
legend(x=0.7, y=1, legend = rownames(data[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "grey", cex=1.2, pt.cex=3)
 
 

 
 
 
 
compare=merge(table1,table2,by="terms")
compare=merge(compare,table3,by="terms")
compare=merge(compare,table4,by="terms")
compare=merge(compare,table5,by="terms")
compare=merge(compare,table6,by="terms")
compare=merge(compare,table7,by="terms")
compare2=compare[,2:8]
compare2=as.matrix(compare2)
compare2=apply(compare2, 2, as.numeric)
row.names(compare2)=compare[,1]
compareWSB=compare2[rowSums(compare2)>0,]


#### Radar plot example

library(fmsb)


compareP=t(compareWSB)
data=as.data.frame(compareP)

fb=data.frame(max=rep(max(data),ncol(data)),min=rep(min(data),ncol(data)))
fb=t(fb)
colnames(fb)=colnames(data)
data=rbind(fb,data)

colors_border=c( rgb(0.5,0.5,0.5,0.7),rgb(1,0.54,0,0.7),rgb(1,0.54,0,0.7), rgb(1,0.54,0,0.7), rgb(0,0.4,0,0.7),rgb(0,0.4,0,0.7),rgb(0,0.4,0,0.7)  )
## colors_in=c( rgb(0.5,0.5,0.5,0.4),rgb(1,0.54,0,0.4), rgb(1,0.54,0,0.4), rgb(0,0.4,0,0.4),rgb(0,0.4,0,0.4),rgb(0,0.4,0,0.4)  )
colors_in=c( rgb(0.5,0.5,0.5,0),rgb(1,0.54,0,0),rgb(1,0.54,0,0), rgb(1,0.54,0,0), rgb(0,0.4,0,0),rgb(0,0.4,0,0),rgb(0,0.4,0,0)  )
radarchart( data  , axistype=1 , 
    #custom polygon
    pcol=colors_border , pfcol=colors_in , plwd=2 , plty=c(1,1:3,1:3),
    #custom the grid
    cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,20,5), cglwd=0.8,
    #custom labels
    vlcex=0.8 
    )
	
legend(x=1, y=1, legend = rownames(data[-c(1,2),]), col=colors_border, text.col = "black", cex=1.2, lwd=2, lty=c(1,1:3,1:3))
 
 




 
 
 
 
 
 
 
 







