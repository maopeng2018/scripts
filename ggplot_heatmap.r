setwd("D:\\Projects\\creA_project\\analysis\\visualize")

library(ggplot2)
library(reshape2)
library(plyr)
library("ggdendro")

data <- read.table("all_cazyme_fold2_p001.txt",sep="\t",head=T)

##nba=data[1:50,1:14]
nba=data[,1:14]
values=nba[,2:14]
nba[,2:14]=log2(nba[,2:14])
rownames(values)=nba$Name
dendro <- as.dendrogram(hclust(d = dist(x = values)))
dendro.plot <- ggdendrogram(data = dendro, rotate = TRUE)
print(dendro.plot)


# Extract the order of the tips in the dendrogram
orders <- order.dendrogram(dendro)
# Order the levels according to their position in the cluster

nba$Name = factor(nba$Name,nba$Name[orders])


nba.m <- melt(nba)
## nba.m <- ddply(nba.m, .(variable), transform, rescale = scale(value))


heatmap <- (ggplot(nba.m, aes(variable, Name)) + geom_tile(aes(fill = value), colour = "white") 
      + scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
## 	  + geom_text(aes(label = round(nba.m$value, 2)),size=1)
	  +theme(axis.text.x = element_text(colour="grey20",size=2,angle=90,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=2,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=2,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=2,angle=90,hjust=.5,vjust=.5,face="plain"))
  )

print(heatmap)





heatmap <- (ggplot(nba.m, aes(variable, Name)) + geom_tile(aes(fill = value), colour = "white") 
      + scale_fill_gradient(low = "green", high = "red")
	  + geom_text(aes(label = round(nba.m$value, 2)),size=1)
	  +theme(axis.text.x = element_text(colour="grey20",size=2,angle=90,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=2,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=2,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=2,angle=90,hjust=.5,vjust=.5,face="plain"))
  )

print(heatmap)


heatmap <- (ggplot(nba.m, aes(variable, Name)) + geom_tile(aes(fill = value), colour = "white") 
      + scale_fill_gradient2()
	  + geom_text(aes(label = round(nba.m$value, 2)),size=1)
	  +theme(axis.text.x = element_text(colour="grey20",size=2,angle=90,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=2,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=2,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=2,angle=90,hjust=.5,vjust=.5,face="plain"))
  )

print(heatmap)

	