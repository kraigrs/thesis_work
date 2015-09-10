data1 <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/imprinting/chr3R/zhr_z30.chr3R_k50.imprinting.txt",header=TRUE,sep="\t");

plot((data1$right_bound+data1$left_bound)/(2*1000000),data1$imprint/50,xlab="Window center (in Mb)",ylab="% imprinted genes",main="Imprinting on chr3R",type="b",pch=19,cex=0.4,ylim=c(0,1),xlim=c(0,28),col=rgb(0,0,0,0.2));

legend("topright",legend="Clusters of 50 genes",fill="black");

##########################################################
# make plot with number of heterozygous sites per window #
##########################################################

par(mar=c(5, 4, 4, 4) + 0.1);

plot((data1$right_bound+data1$left_bound)/(2*1000000),data1$imprint/50,xlab="Window center (in Mb)",ylab="% imprinted genes",main="Imprinting on chr3R",type="b",pch=19,cex=0.4,ylim=c(0,1),xlim=c(0,28),col=rgb(0,0,0,0.2));

par(new=T);

temp1 <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/imprinting/chr3R/chr3R_hetSNPs.bed",sep="\t");
temp2 <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/imprinting/chr3R/chr3R_hetINDELs.bed",sep="\t");
temp3 <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/imprinting/chr3R/chr3R_hetSNPs.bed",sep="\t");
temp4 <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/imprinting/chr3R/chr3R_hetINDELs.bed",sep="\t");

tots <- temp1[,4]+temp2[,4]+temp3[,4]+temp4[,4];

temp <- cbind(temp1,tots);

plot((temp[,2]+temp[,3])/(2*1000000),temp[,5],ylim=c(0,15000),xlim=c(0,28),axes=F,xlab="",ylab="",main="",type="l",col=rgb(0,0,1,0.2),pch=19,cex=0.3);

axis(4, ylim=c(0,15000), col="black",col.axis="black");
mtext("# heterozygous sites (SNPs and INDELs)",side=4,col="black",line=2.5);

legend("topleft",legend="Clusters of 50 genes",fill="black");

##################################
# make plot with -log10(pvalues) #
##################################

gene_pvals_3R <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/null_dist/chr3R/zhr_z30.chr3R_k50.imprinting_pvals.txt",header=TRUE,sep="\t");
gene_pvals_3L <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/null_dist/chr3L/zhr_z30.chr3L_k50.imprinting_pvals.txt",header=TRUE,sep="\t");
gene_pvals_2R <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/null_dist/chr2R/zhr_z30.chr2R_k50.imprinting_pvals.txt",header=TRUE,sep="\t");
gene_pvals_2L <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/null_dist/chr2L/zhr_z30.chr2L_k50.imprinting_pvals.txt",header=TRUE,sep="\t");
gene_pvals_X <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/null_dist/chrX/zhr_z30.chrX_k50.imprinting_pvals.txt",header=TRUE,sep="\t");

gene_pvals <- rbind(gene_pvals_3R, gene_pvals_3L, gene_pvals_2R, gene_pvals_2L, gene_pvals_X);
gene_qvals <- p.adjust(gene_pvals$pval,method="fdr");

temp <- cbind(gene_pvals,gene_qvals);
temp <- subset(temp,temp$chr == "3R");
temp <- temp[order(temp$left_bound),];

par(mar=c(5, 4, 4, 4) + 0.1);

plot((data1$right_bound+data1$left_bound)/(2*1000000),data1$imprint/50,xlab="Window center (in Mb)",ylab="% imprinted genes",main="Imprinting on chr3R",type="b",pch=19,cex=0.4,ylim=c(0,1),xlim=c(0,28),col=rgb(0,0,0,0.2));

par(new=T);

plot((temp[,2]+temp[,3])/(2*1000000),-log10(temp[,5]),ylim=c(0,4.5),xlim=c(0,28),axes=F,xlab="",ylab="",main="",type="l",col=rgb(0,0,1,0.2),pch=19,cex=0.3);

axis(4, ylim=c(0,3), col="black",col.axis="black");
mtext(text=expression(paste(-log[10],"(q-value) for enrichment of imprinting",sep="")),side=4,col="black",line=2.5);

legend("topleft",legend="Clusters of 50 genes",fill="black",bty="n");
legend("topright",legend=c("q-value = 0.05"),lty=2,bty="n");
abline(h=-log10(0.05),lty=2);