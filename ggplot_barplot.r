setwd("D:\\Projects\\creA_project\\analysis\\Target_gene_analysis\\CreA_motif");
data=read.table("KOG_distribution.txt",sep="\t",head=T)

library(ggplot2)

data2 <- data.frame(
    group = rep(c("A.niger genome","Genes with CreA motif","Genes with CreA motif and derepression"), each=25),
    x = rep(data$terms, 3),
    y = 100*(c(as.vector(data[,2]),as.vector(data[,3]),as.vector(data[,4])))
)



  
ggplot(data2, aes(x=reorder(x,1:75), y=y, fill=group)) + 
  geom_bar(stat="identity", position="dodge") +
  theme(axis.text.x = element_text(face="bold", size=10, angle=45))+
  scale_fill_grey(start = 0.8, end = 0.2)  ##  scale_fill_brewer(palette="Pastel1")



