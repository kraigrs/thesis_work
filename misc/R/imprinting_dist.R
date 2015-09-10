data2 <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/imprinting/chr3R/zhr_z30.chr3R_wind500000_step10000.imprinting.txt",header=TRUE,sep="\t");

#############################################
# make plot with number of genes per window #
#############################################

par(mar=c(5, 4, 4, 4) + 0.1);

plot((data2$right_bound+data2$left_bound)/(2*1000000),data2$imprinted/data2$total,xlab="Window center (in Mb)",ylab="% imprinted genes",main="Imprinting on chr3R",type="p",pch=19,cex=0.4,ylim=c(0,1),xlim=c(0,28),col=rgb(0,0,0,0.2));

par(new=T);

temp <- data2[order(data2$left_bound),];

plot((temp$right_bound+temp$left_bound)/(2*1000000),temp$total,ylim=c(0,100),xlim=c(0,28),axes=F,xlab="",ylab="",main="",type="l",col=rgb(0,0,1,0.2),pch=19,cex=0.3);

axis(4, ylim=c(0,100), col="black",col.axis="black");
mtext("# genes in window",side=4,col="black",line=2.5);

legend("topright",legend="500Kb windows with 10Kb offset",fill="black");

##########################################################
# make plot with number of heterozygous sites per window #
##########################################################

par(mar=c(5, 4, 4, 4) + 0.1);

plot((data2$right_bound+data2$left_bound)/(2*1000000),data2$imprinted/data2$total,xlab="Window center (in Mb)",ylab="% imprinted genes",main="Imprinting on chr3R",type="p",pch=19,cex=0.4,ylim=c(0,1),xlim=c(0,28),col=rgb(0,0,0,0.2));

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

legend("topleft",legend="500Kb windows with 10Kb offset",fill="black");

##################################
# make plot with -log10(pvalues) #
##################################

dist_pvals_3R <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/null_dist/chr3R/zhr_z30.chr3R_wind500000_step10000.imprinting_pvals.txt",header=TRUE,sep="\t");
dist_pvals_3L <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/null_dist/chr3L/zhr_z30.chr3L_wind500000_step10000.imprinting_pvals.txt",header=TRUE,sep="\t");
dist_pvals_2R <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/null_dist/chr2R/zhr_z30.chr2R_wind500000_step10000.imprinting_pvals.txt",header=TRUE,sep="\t");
dist_pvals_2L <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/null_dist/chr2L/zhr_z30.chr2L_wind500000_step10000.imprinting_pvals.txt",header=TRUE,sep="\t");
dist_pvals_X <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/null_dist/chrX/zhr_z30.chrX_wind500000_step10000.imprinting_pvals.txt",header=TRUE,sep="\t");
dist_pvals_4 <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/null_dist/chr4/zhr_z30.chr4_wind500000_step10000.imprinting_pvals.txt",header=TRUE,sep="\t");

dist_pvals <- rbind(dist_pvals_3R, dist_pvals_3L, dist_pvals_2R, dist_pvals_2L, dist_pvals_X, dist_pvals_4);
dist_qvals <- p.adjust(dist_pvals$pval,method="fdr");
dist_sig <- cbind(dist_pvals,dist_qvals);

temp <- cbind(dist_pvals,dist_qvals);
temp <- subset(temp,temp$chr == "3R");
temp <- temp[order(temp$left_bound),];

par(mar=c(5, 4, 4, 4) + 0.1);

plot((data2$right_bound+data2$left_bound)/(2*1000000),data2$imprinted/data2$total,xlab="Window center (in Mb)",ylab="% imprinted genes",main="Imprinting on chr3R",type="p",pch=19,cex=0.4,ylim=c(0,1),xlim=c(0,28),col=rgb(0,0,0,0.2));

par(new=T);

plot((temp[,2]+temp[,3])/(2*1000000),-log10(temp[,5]),ylim=c(0,4.5),xlim=c(0,28),axes=F,xlab="",ylab="",main="",type="l",col=rgb(0,0,1,0.2),pch=19,cex=0.3);

axis(4, ylim=c(0,3), col="black",col.axis="black");
mtext(text=expression(paste(-log[10],"(q-value) for enrichment of imprinting",sep="")),side=4,col="black",line=2.5);

legend("topleft",legend="500Kb windows with 10Kb offset",fill="black",bty="n");
legend("topright",legend=c("q-value = 0.05"),lty=2,bty="n");
abline(h=-log10(0.05),lty=2);

###############################################
# plot with -log10(pvalues) for 6 chromosomes #
###############################################

pdf(file="/Users/kraigrs/Desktop/color_imprinting_plot.pdf",width=10,height=7);

chrX_dist <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/imprinting/chrX/zhr_z30.chrX_wind500000_step10000.imprinting.txt",header=TRUE,sep="\t");
chr2R_dist <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/imprinting/chr2R/zhr_z30.chr2R_wind500000_step10000.imprinting.txt",header=TRUE,sep="\t");
chr2L_dist <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/imprinting/chr2L/zhr_z30.chr2L_wind500000_step10000.imprinting.txt",header=TRUE,sep="\t");
chr3R_dist <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/imprinting/chr3R/zhr_z30.chr3R_wind500000_step10000.imprinting.txt",header=TRUE,sep="\t");
chr3L_dist <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/imprinting/chr3L/zhr_z30.chr3L_wind500000_step10000.imprinting.txt",header=TRUE,sep="\t");
chr4_dist <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/imprinting/chr4/zhr_z30.chr4_wind500000_step10000.imprinting.txt",header=TRUE,sep="\t");

chr <- rep("X",nrow(chrX_dist));
chrX_dist <- cbind(chrX_dist,chr);

chr <- rep("3R",nrow(chr3R_dist));
chr3R_dist <- cbind(chr3R_dist,chr);

chr <- rep("3L",nrow(chr3L_dist));
chr3L_dist <- cbind(chr3L_dist,chr);

chr <- rep("2R",nrow(chr2R_dist));
chr2R_dist <- cbind(chr2R_dist,chr);

chr <- rep("2L",nrow(chr2L_dist));
chr2L_dist <- cbind(chr2L_dist,chr);

chr <- rep("4LR",nrow(chr4_dist));
chr4_dist <- cbind(chr4_dist,chr);

chr_dist <- rbind(chrX_dist,chr2L_dist,chr2R_dist,chr3L_dist,chr3R_dist,chr4_dist);

dist_pvals_3R <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/null_dist/chr3R/zhr_z30.chr3R_wind500000_step10000.imprinting_pvals.txt",header=TRUE,sep="\t");
dist_pvals_3L <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/null_dist/chr3L/zhr_z30.chr3L_wind500000_step10000.imprinting_pvals.txt",header=TRUE,sep="\t");
dist_pvals_2R <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/null_dist/chr2R/zhr_z30.chr2R_wind500000_step10000.imprinting_pvals.txt",header=TRUE,sep="\t");
dist_pvals_2L <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/null_dist/chr2L/zhr_z30.chr2L_wind500000_step10000.imprinting_pvals.txt",header=TRUE,sep="\t");
dist_pvals_X <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/null_dist/chrX/zhr_z30.chrX_wind500000_step10000.imprinting_pvals.txt",header=TRUE,sep="\t");
dist_pvals_4 <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/null_dist/chr4/zhr_z30.chr4_wind500000_step10000.imprinting_pvals.txt",header=TRUE,sep="\t");

dist_pvals <- rbind(dist_pvals_3R, dist_pvals_3L, dist_pvals_2R, dist_pvals_2L, dist_pvals_X, dist_pvals_4);
dist_qvals <- p.adjust(dist_pvals$pval,method="fdr");
dist_sig <- cbind(dist_pvals,dist_qvals);

dist_data <- merge(chr_dist,dist_sig);

library(lattice);

#xyplot( imprinted/total + total ~ (right_bound+left_bound)/(2*1000000) | chr, dist_data);

###################
#### Manhattan plot
###################

##### imprinted / total genes ######
par(mar=c(5, 4, 4, 4) + 0.1);

plot((chrX_dist$right_bound+chrX_dist$left_bound)/(2*1000000),chrX_dist$imprinted/chrX_dist$total*100,xlab="",ylab="",main="",type="p",pch=19,cex=0.4,ylim=c(0,100),xlim=c(0,120),col="red",xaxt="n",yaxt="n");

par(new=T);
plot((chr2L_dist$right_bound+chr2L_dist$left_bound)/(2*1000000)+max((chrX_dist$right_bound+chrX_dist$left_bound)/(2*1000000))+0.000001,chr2L_dist$imprinted/chr2L_dist$total*100,xlab="",ylab="",main="",type="p",pch=19,cex=0.4,ylim=c(0,100),xlim=c(0,120),col="blue",xaxt="n",yaxt="n");

par(new=T);
plot((chr2R_dist$right_bound+chr2R_dist$left_bound)/(2*1000000)+max((chrX_dist$right_bound+chrX_dist$left_bound)/(2*1000000))+max((chr2L_dist$right_bound+chr2L_dist$left_bound)/(2*1000000))+0.000001,chr2R_dist$imprinted/chr2R_dist$total*100,xlab="",ylab="",main="",type="p",pch=19,cex=0.4,ylim=c(0,100),xlim=c(0,120),col="forestgreen",xaxt="n",yaxt="n");

par(new=T);
plot((chr3L_dist$right_bound+chr3L_dist$left_bound)/(2*1000000)+max((chrX_dist$right_bound+chrX_dist$left_bound)/(2*1000000))+max((chr2L_dist$right_bound+chr2L_dist$left_bound)/(2*1000000))+max((chr2R_dist$right_bound+chr2R_dist$left_bound)/(2*1000000))+0.000001,chr3L_dist$imprinted/chr3L_dist$total*100,xlab="",ylab="",main="",type="p",pch=19,cex=0.4,ylim=c(0,100),xlim=c(0,120),col="orange",xaxt="n",yaxt="n");

par(new=T);
plot((chr3R_dist$right_bound+chr3R_dist$left_bound)/(2*1000000)+max((chrX_dist$right_bound+chrX_dist$left_bound)/(2*1000000))+max((chr2L_dist$right_bound+chr2L_dist$left_bound)/(2*1000000))+max((chr2R_dist$right_bound+chr2R_dist$left_bound)/(2*1000000))+max((chr3L_dist$right_bound+chr3L_dist$left_bound)/(2*1000000))+0.000001,chr3R_dist$imprinted/chr3R_dist$total*100,xlab="",ylab="",main="",type="p",pch=19,cex=0.4,ylim=c(0,100),xlim=c(0,120),col="purple",xaxt="n",yaxt="n");

par(new=T);
plot((chr4_dist$right_bound+chr4_dist$left_bound)/(2*1000000)+max((chrX_dist$right_bound+chrX_dist$left_bound)/(2*1000000))+max((chr2L_dist$right_bound+chr2L_dist$left_bound)/(2*1000000))+max((chr2R_dist$right_bound+chr2R_dist$left_bound)/(2*1000000))+max((chr3L_dist$right_bound+chr3L_dist$left_bound)/(2*1000000))+max((chr3R_dist$right_bound+chr3R_dist$left_bound)/(2*1000000))+0.000001,chr4_dist$imprinted/chr4_dist$total*100,xlab="",main="Significant imprinting in clusters along genome",type="p",pch=19,cex=0.4,ylim=c(0,100),xlim=c(0,120),col="yellow",xaxt="n",ylab="% imprinted genes");

##### q-values#######

temp1 <- subset(dist_data,dist_data$chr == "X");
temp1 <- temp1[order(temp1$left_bound),];
par(new=T);
plot((temp1$right_bound+temp1$left_bound)/(2*1000000),-log10(temp1$dist_qvals),xlab="",ylab="",main="",type="l",pch=19,cex=0.4,ylim=c(0,2),xlim=c(0,120),col=rgb(0,0,0,0.5),axes=F,xaxt="n");

temp2 <- subset(dist_data,dist_data$chr == "2L");
temp2 <- temp2[order(temp2$left_bound),];
par(new=T);
plot((temp2$right_bound+temp2$left_bound)/(2*1000000)+max((temp1$right_bound+temp1$left_bound)/(2*1000000))+0.000001,-log10(temp2$dist_qvals),xlab="",ylab="",main="",type="l",pch=19,cex=0.4,ylim=c(0,2),xlim=c(0,120),col=rgb(0,0,0,0.5),axes=F,xaxt="n");

temp3 <- subset(dist_data,dist_data$chr == "2R");
temp3 <- temp3[order(temp3$left_bound),];
par(new=T);
plot((temp3$right_bound+temp3$left_bound)/(2*1000000)+max((temp2$right_bound+temp2$left_bound)/(2*1000000))+max((temp1$right_bound+temp1$left_bound)/(2*1000000))+0.000001,-log10(temp3$dist_qvals),xlab="",ylab="",main="",type="l",pch=19,cex=0.4,ylim=c(0,2),xlim=c(0,120),col=rgb(0,0,0,0.5),axes=F,xaxt="n");

temp4 <- subset(dist_data,dist_data$chr == "3L");
temp4 <- temp4[order(temp4$left_bound),];
par(new=T);
plot((temp4$right_bound+temp4$left_bound)/(2*1000000)+max((temp3$right_bound+temp3$left_bound)/(2*1000000))+max((temp2$right_bound+temp2$left_bound)/(2*1000000))+max((temp1$right_bound+temp1$left_bound)/(2*1000000))+0.000001,-log10(temp4$dist_qvals),xlab="",ylab="",main="",type="l",pch=19,cex=0.4,ylim=c(0,2),xlim=c(0,120),col=rgb(0,0,0,0.5),axes=F,xaxt="n");

temp5 <- subset(dist_data,dist_data$chr == "3R");
temp5 <- temp5[order(temp5$left_bound),];
par(new=T);
plot((temp5$right_bound+temp5$left_bound)/(2*1000000)+max((temp4$right_bound+temp4$left_bound)/(2*1000000))+max((temp3$right_bound+temp3$left_bound)/(2*1000000))+max((temp2$right_bound+temp2$left_bound)/(2*1000000))+max((temp1$right_bound+temp1$left_bound)/(2*1000000))+0.000001,-log10(temp5$dist_qvals),xlab="",ylab="",main="",type="l",pch=19,cex=0.4,ylim=c(0,2),xlim=c(0,120),col=rgb(0,0,0,0.5),axes=F,xaxt="n");

temp6 <- subset(dist_data,dist_data$chr == "4LR");
temp6 <- temp6[order(temp6$left_bound),];
par(new=T);
plot((temp6$right_bound+temp6$left_bound)/(2*1000000)+max((temp5$right_bound+temp5$left_bound)/(2*1000000))+max((temp4$right_bound+temp4$left_bound)/(2*1000000))+max((temp3$right_bound+temp3$left_bound)/(2*1000000))+max((temp2$right_bound+temp2$left_bound)/(2*1000000))+max((temp1$right_bound+temp1$left_bound)/(2*1000000))+0.000001,-log10(temp6$dist_qvals),xlab="",ylab="",main="",type="l",pch=19,cex=0.4,ylim=c(0,2),xlim=c(0,120),col=rgb(0,0,0,0.5),axes=F,xaxt="n");

axis(4, ylim=c(0,3), col="black",col.axis="black");
mtext(text=expression(paste(-log[10],"(q-value) for enrichment of imprinting",sep="")),side=4,col="black",line=2.5);

#legend("topleft",legend=c("X","2L","2R","3L","3R","4"),fill=c("red","blue","forestgreen","orange","purple","yellow"),bty="n");
#legend("top",legend=c("q-value = 0.05"),lty=2,bty="n");
abline(h=-log10(0.05),lty=2);

dev.off();



