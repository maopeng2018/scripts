setwd("D:\\test\\heatmap\\Ronald")
getwd();

library("gplots")

## cazyFPKM=read.table("select_transporter_candidates_expression.txt",sep="\t",head=T)
data=read.table("correlation file.txt",sep="\t",head=T)

fpkm=data[,2:ncol(data)]
rownames(fpkm)=data[,1]
values2=fpkm

values2=as.matrix(values2)
 
windows()
col2=topo.colors(75)[75:1]
col3=col2[c(1,15,25,30,36,42,48,50:75)]
dist.pear <- function(x) as.dist((1-cor(t(x),method="pearson"))/2)
hclust.ave <- function(x) hclust(x, method="complete")
heatmap.2(values2, margin=c(8,15), cexCol=0.8,cexRow=0.6, density.info="none", col=col3, trace="none", distfun=dist.pear, hclustfun=hclust.ave)

		 

windows()
dist.pear <- function(x) as.dist((1-cor(t(x),method="spearman")))
hclust.ave <- function(x) hclust(x, method="complete")
heatmap.2(values2, density.info="none", col=topo.colors(75), trace="none", distfun=dist.pear, hclustfun=hclust.ave)
		 
windows()
dist.pear <- function(x) as.dist((1-cor(t(x),method="pearson"))/2)
hclust.ave <- function(x) hclust(x, method="complete")
heatmap.2(values2, margin=c(8,15), cexCol=0.8,cexRow=0.6, density.info="none", col=topo.colors(75), trace="none", distfun=dist.pear, hclustfun=hclust.ave)




windows()
col2=topo.colors(20)[20:1]
dist.pear <- function(x) as.dist((1-cor(t(x),method="pearson"))/2)
hclust.ave <- function(x) hclust(x, method="complete")
heatmap.2(values2, margin=c(8,15), cexCol=0.8,cexRow=0.6, density.info="none", col=col2, trace="none", distfun=dist.pear, hclustfun=hclust.ave)


windows()
library("RColorBrewer")
coul = brewer.pal(9, "YlGnBu")
dist.pear <- function(x) as.dist((1-cor(t(x),method="pearson"))/2)
hclust.ave <- function(x) hclust(x, method="complete")
heatmap.2(values2, margin=c(8,15), cexCol=0.8,cexRow=0.6, density.info="none", col=coul, trace="none", distfun=dist.pear, hclustfun=hclust.ave)



windows()
heatmap.2(values2,
           margin=c(8,15),
		     col=topo.colors(75),
            ## Colv=NA,
			 trace="none",
			 cexCol=0.8,
			 cexRow=0.6,
			 key=T,
			 density.info="none",
			 keysize=0.75,
			 distfun = function(x) dist(x,method = 'euclidean'),    ### "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski" 
			hclustfun = function(x) hclust(x,method = 'complete'),   ###  other methods including "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).'centroid',‘complete’
           ##   ylab="Plant biomass degeradation Cazy genes",			
          ) 

		  





		  
		 
heatmap.2(values2,
            margin=c(5,15),
col=bluered(75),
scale="none", 
 key=TRUE, symkey=FALSE, density.info="none", 
  trace="none",
 cexCol=1,
 cexRow=0.3,
      ##ylab="Plant biomass degeradation Cazy genes",
         ) 


		 
library(RColorBrewer)
coul = brewer.pal(4, "BuPu")
heatmap.2(values2,
            margin=c(5,15),
col= colorRampPalette(coul)(25),
scale="none", 
 key=TRUE, symkey=FALSE, density.info="none", 
  trace="none",
 cexCol=1,
 cexRow=0.3,
      ##ylab="Plant biomass degeradation Cazy genes",
         ) 
		 

library(RColorBrewer)
coul = brewer.pal(11, "RdBu")
## display.brewer.pal(n = 11, name = 'RdBu')
heatmap.2(values2,
            margin=c(5,15),
col= colorRampPalette(coul)(11),
scale="none", 
 key=TRUE, symkey=FALSE, density.info="none", 
  trace="none",
 cexCol=1,
 cexRow=0.3,
      ##ylab="Plant biomass degeradation Cazy genes",
         ) 
		 
library(RColorBrewer)
coul = brewer.pal(11, "RdBu")
## display.brewer.pal(n = 11, name = 'RdBu')
heatmap.2(values2,
            margin=c(5,15),
col= colorRampPalette(coul)(11),
scale="none", 
 key=TRUE, symkey=FALSE, density.info="none", 
  trace="none",
 cexCol=1,
 cexRow=0.3,
      ##ylab="Plant biomass degeradation Cazy genes",
         ) 
		 