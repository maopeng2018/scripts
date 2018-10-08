setwd("D:\\Projects\\Evy_mycelicphlicia\\Daniel_data\\statis_analysis\\target_analysis")
getwd();


metabolic=read.table("PPP_PCP_genes.txt",sep="\t",head=T)

RPKM=read.table("expression_RPKM.txt",sep="\t",head=T)

metab=merge(metabolic,RPKM,by.x="JGI_id",by.y="GeneID")

write.table(metab,file="All_metabolizeGene_expression_profiles.txt",sep="\t",col.names=NA)


mean_gluW=apply(metab[,6:7],1,mean)
mean_gluM=apply(metab[,8:9],1,mean)
mean_araW=apply(metab[,10:11],1,mean)
mean_araM=apply(metab[,12:13],1,mean)
mean_xylW=apply(metab[,14:15],1,mean)
mean_xylM=apply(metab[,16:17],1,mean)

exp2=cbind(mean_gluW,mean_gluM,mean_araW,mean_araM,mean_xylW,mean_xylM);
exp=log2(exp2+1)


library(RColorBrewer)
colour= colorRampPalette(brewer.pal(5, "Blues"))(25)


genes=(dim(exp))[1]

mat1 <- matrix(exp[1,], nrow = 6, ncol = 1)
mat2 <- matrix(exp[2,], nrow = 6, ncol = 1)
mat3 <- matrix(exp[3,], nrow = 6, ncol = 1)
mat4 <- matrix(exp[4,], nrow = 6, ncol = 1)
mat5 <- matrix(exp[5,], nrow = 6, ncol = 1)
mat6 <- matrix(exp[6,], nrow = 6, ncol = 1)
mat7 <- matrix(exp[7,], nrow = 6, ncol = 1)
mat8 <- matrix(exp[8,], nrow = 6, ncol = 1)
mat9 <- matrix(exp[9,], nrow = 6, ncol = 1)
mat10 <- matrix(exp[10,], nrow = 6, ncol = 1)
mat11 <- matrix(exp[11,], nrow = 6, ncol = 1)
mat12 <- matrix(exp[12,], nrow = 6, ncol = 1)
mat13 <- matrix(exp[13,], nrow = 6, ncol = 1)
mat14 <- matrix(exp[14,], nrow = 6, ncol = 1)
mat15 <- matrix(exp[15,], nrow = 6, ncol = 1)
mat16 <- matrix(exp[16,], nrow = 6, ncol = 1)
mat17 <- matrix(exp[17,], nrow = 6, ncol = 1)
mat18 <- matrix(exp[18,], nrow = 6, ncol = 1)
mat19 <- matrix(exp[19,], nrow = 6, ncol = 1)
mat20 <- matrix(exp[20,], nrow = 6, ncol = 1)


zz <- list(mat1,mat2,mat3,mat4,mat5,mat6,mat7,mat8,mat9,mat10,mat11,mat12,mat13,mat14,mat15,mat16,mat17,mat18,mat19,mat20)

par(mfrow=c(5,1))
z=unlist(zz) 
a=round(min(z));  b=round(max(z));   c=b-a+1;
m=matrix(seq(a,b,by=0.5),c*2-1,1)
image(m,col=colour,zlim=range(a,b)) #### here 0=mi


par(mfrow=c(5,1))  ### show multiple images in a figure
for(i in 1:5){
image(zz[[i]],col=colour,zlim=range(a,b),axe=F)       #### zlim=is crucial
abline(v=c(0.1,0.3,0.5,0.7,0.9), col="black", lwd=1)
}

par(mfrow=c(5,1))  ### show multiple images in a figure
for(i in 6:10){
image(zz[[i]],col=colour,zlim=range(a,b),axe=F)       #### zlim=is crucial
abline(v=c(0.1,0.3,0.5,0.7,0.9), col="black", lwd=1)
}

par(mfrow=c(5,1))  ### show multiple images in a figure
for(i in 11:15){
image(zz[[i]],col=colour,zlim=range(a,b),axe=F)       #### zlim=is crucial
abline(v=c(0.1,0.3,0.5,0.7,0.9), col="black", lwd=1)
}

par(mfrow=c(5,1))  ### show multiple images in a figure
for(i in 16:20){
image(zz[[i]],col=colour,zlim=range(a,b),axe=F)       #### zlim=is crucial
abline(v=c(0.1,0.3,0.5,0.7,0.9), col="black", lwd=1)
}

par(mfrow=c(5,1))  ### show multiple images in a figure
for(i in 14:18){
image(zz[[i]],col=colour,zlim=range(a,b),axe=F)       #### zlim=is crucial
abline(v=c(0.1,0.3,0.5,0.7,0.9), col="black", lwd=1)
}


