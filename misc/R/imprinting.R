data1 <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/imprinting/chr3R/zhr_z30.chr3R_k50.imprinting.txt",header=TRUE,sep="\t");

data2 <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/imprinting/chr3R/zhr_z30.chr3R_wind500000_step10000.imprinting.txt",header=TRUE,sep="\t");

#############################################
# make plot with number of genes per window #
#############################################

par(mar=c(5, 4, 4, 4) + 0.1);

plot((data1$right_bound+data1$left_bound)/(2*1000000),data1$imprint/50,xlab="Window center (in Mb)",ylab="% imprinted genes",main="Imprinting on chr3R",type="b",pch=19,cex=0.4,ylim=c(0,1),xlim=c(0,28),col=rgb(0,0,0,0.2));

par(new=T);

plot((data2$right_bound+data2$left_bound)/(2*1000000),data2$imprint/data2$total,xlab="Window center (in Mb)",ylab="% imprinted genes",main="Imprinting on chr3R",type="p",pch=19,cex=0.4,ylim=c(0,1),xlim=c(0,28),col=rgb(1,0,0,0.2));

par(new=T);

temp <- data2[order(data2$left_bound),];

plot((temp$right_bound+temp$left_bound)/(2*1000000),temp$total,ylim=c(0,100),xlim=c(0,28),axes=F,xlab="",ylab="",main="",type="l",col=rgb(0,0,1,0.2),pch=19,cex=0.3);

axis(4, ylim=c(0,100), col="black",col.axis="black");
mtext("# genes in window",side=4,col="black",line=2.5);

legend("topright",legend=c("Clusters of 50 genes","500Kb windows with 10Kb offset"),fill=c("black","red"));
#############################################

##########################################################
# make plot with number of heterozygous sites per window #
##########################################################

par(mar=c(5, 4, 4, 4) + 0.1);

plot((data1$right_bound+data1$left_bound)/(2*1000000),data1$imprint/50,xlab="Window center (in Mb)",ylab="% imprinted genes",main="Imprinting on chr3R",type="b",pch=19,cex=0.4,ylim=c(0,1),xlim=c(0,28),col=rgb(0,0,0,0.2));

par(new=T);

plot((data2$right_bound+data2$left_bound)/(2*1000000),data2$imprint/data2$total,xlab="Window center (in Mb)",ylab="% imprinted genes",main="Imprinting on chr3R",type="p",pch=19,cex=0.4,ylim=c(0,1),xlim=c(0,28),col=rgb(1,0,0,0.2));

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

legend("topleft",legend=c("Clusters of 50 genes","500Kb windows with 10Kb offset"),fill=c("black","red"));

##########################################################

##################################
# make plot with -log10(pvalues) #
##################################

dist_pvals_3R <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/null_dist/chr3R/zhr_z30.chr3R_wind500000_step10000.imprinting_pvals.txt",header=TRUE,sep="\t");
dist_pvals_3L <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/null_dist/chr3L/zhr_z30.chr3L_wind500000_step10000.imprinting_pvals.txt",header=TRUE,sep="\t");
dist_pvals_2R <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/null_dist/chr2R/zhr_z30.chr2R_wind500000_step10000.imprinting_pvals.txt",header=TRUE,sep="\t");
dist_pvals_2L <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/null_dist/chr2L/zhr_z30.chr2L_wind500000_step10000.imprinting_pvals.txt",header=TRUE,sep="\t");
dist_pvals_X <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/null_dist/chrX/zhr_z30.chrX_wind500000_step10000.imprinting_pvals.txt",header=TRUE,sep="\t");

dist_pvals <- rbind(dist_pvals_3R, dist_pvals_3L, dist_pvals_2R, dist_pvals_2L, dist_pvals_X);

x_dist <- runif(nrow(dist_pvals));
dist_qvals <- p.adjust(dist_pvals$pval,method="fdr");

data1 <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/imprinting/chr3R/zhr_z30.chr3R_k50.imprinting.txt",header=TRUE,sep="\t");
data2 <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/imprinting/chr3R/zhr_z30.chr3R_wind500000_step10000.imprinting.txt",header=TRUE,sep="\t");

temp <- cbind(dist_pvals,dist_qvals);
temp <- subset(temp,temp$chr == "3R");
temp <- temp[order(temp$left_bound),];

par(mar=c(5, 4, 4, 4) + 0.1);

plot((data1$right_bound+data1$left_bound)/(2*1000000),data1$imprint/50,xlab="",ylab="",main="",type="b",pch=19,cex=0.4,ylim=c(0,1),xlim=c(0,28),col=rgb(0,0,0,0.2));

par(new=T);

plot((data2$right_bound+data2$left_bound)/(2*1000000),data2$imprint/data2$total,xlab="Window center (in Mb)",ylab="% imprinted genes",main="Imprinting on chr3R",type="p",pch=19,cex=0.4,ylim=c(0,1),xlim=c(0,28),col=rgb(1,0,0,0.2));

par(new=T);

plot((temp[,2]+temp[,3])/(2*1000000),-log10(temp[,5]),ylim=c(0,4.5),xlim=c(0,28),axes=F,xlab="",ylab="",main="",type="l",col=rgb(0,0,1,0.2),pch=19,cex=0.3);

axis(4, ylim=c(0,3), col="black",col.axis="black");
mtext(text=expression(paste(-log[10],"(q-value) for enrichment of imprinting",sep="")),side=4,col="black",line=2.5);

legend("topleft",legend=c("Clusters of 50 genes","500Kb windows with 10Kb offset"),fill=c("black","red"),bty="n");
legend("topright",legend=c("q-value = 0.05"),lty=2,bty="n");
abline(h=-log10(0.05),lty=2);
##########################################################







