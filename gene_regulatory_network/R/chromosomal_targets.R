########## expression data ##########

chrom_2L <- read.table("/Users/kraigrs/Wittkopp/gene_regulatory_network/Dmel_chromosomal_targets_2L.txt",header=TRUE,sep="\t");
chrom_2L <- chrom_2L[which(chrom_2L$chromosome != 4),];
chrom_2L$chromosome <- factor(chrom_2L$chromosome,levels=c("X","2L","2R","3L","3R"));

chrom_2R <- read.table("/Users/kraigrs/Wittkopp/gene_regulatory_network/Dmel_chromosomal_targets_2R.txt",header=TRUE,sep="\t");
chrom_2R <- chrom_2R[which(chrom_2R$chromosome != 4),];
chrom_2R$chromosome <- factor(chrom_2R$chromosome,levels=c("X","2L","2R","3L","3R"));

chrom_3L <- read.table("/Users/kraigrs/Wittkopp/gene_regulatory_network/Dmel_chromosomal_targets_3L.txt",header=TRUE,sep="\t");
chrom_3L <- chrom_3L[which(chrom_3L$chromosome != 4),];
chrom_3L$chromosome <- factor(chrom_3L$chromosome,levels=c("X","2L","2R","3L","3R"));

chrom_3R <- read.table("/Users/kraigrs/Wittkopp/gene_regulatory_network/Dmel_chromosomal_targets_3R.txt",header=TRUE,sep="\t");
chrom_3R <- chrom_3R[which(chrom_3R$chromosome != 4),];
chrom_3R$chromosome <- factor(chrom_3R$chromosome,levels=c("X","2L","2R","3L","3R"));

chrom_X <- read.table("/Users/kraigrs/Wittkopp/gene_regulatory_network/Dmel_chromosomal_targets_X.txt",header=TRUE,sep="\t");
chrom_X <- chrom_X[which(chrom_X$chromosome != 4),];
chrom_X$chromosome <- factor(chrom_X$chromosome,levels=c("X","2L","2R","3L","3R"));

pdf(file="/Users/kraigrs/Desktop/Dmel_chromosomal_targets.pdf",height=8.5,width=11);

par(mfrow=c(2,3),oma = c(0,0,3,0));

boxplot(chr2L_regs/regulators ~ chromosome,data=chrom_2L,varwidth=TRUE,ylim=c(0,1),xlab="",ylab="proportion regulators",main="2L");
boxplot(chr3L_regs/regulators ~ chromosome,data=chrom_3L,varwidth=TRUE,ylim=c(0,1),xlab="",ylab="proportion regulators",main="3L");
boxplot(chrX_regs/regulators ~ chromosome,data=chrom_X,varwidth=TRUE,ylim=c(0,1),xlab="chromosome",ylab="proportion regulators",main="X");
boxplot(chr2R_regs/regulators ~ chromosome,data=chrom_2R,varwidth=TRUE,ylim=c(0,1),xlab="chromosome",ylab="proportion regulators",main="2R");
boxplot(chr3R_regs/regulators ~ chromosome,data=chrom_3R,varwidth=TRUE,ylim=c(0,1),xlab="chromosome",ylab="",main="3R");

mtext("expression data", outer = TRUE, cex = 1.5);

dev.off();

########## Kellis data ##########

chrom_2L <- read.table("/Users/kraigrs/Wittkopp/gene_regulatory_network/Kellis_chromosomal_targets_2L.txt",header=TRUE,sep="\t");
chrom_2L <- chrom_2L[which(chrom_2L$chromosome != 4),];
chrom_2L$chromosome <- factor(chrom_2L$chromosome,levels=c("X","2L","2R","3L","3R"));

chrom_2R <- read.table("/Users/kraigrs/Wittkopp/gene_regulatory_network/Kellis_chromosomal_targets_2R.txt",header=TRUE,sep="\t");
chrom_2R <- chrom_2R[which(chrom_2R$chromosome != 4),];
chrom_2R$chromosome <- factor(chrom_2R$chromosome,levels=c("X","2L","2R","3L","3R"));

chrom_3L <- read.table("/Users/kraigrs/Wittkopp/gene_regulatory_network/Kellis_chromosomal_targets_3L.txt",header=TRUE,sep="\t");
chrom_3L <- chrom_3L[which(chrom_3L$chromosome != 4),];
chrom_3L$chromosome <- factor(chrom_3L$chromosome,levels=c("X","2L","2R","3L","3R"));

chrom_3R <- read.table("/Users/kraigrs/Wittkopp/gene_regulatory_network/Kellis_chromosomal_targets_3R.txt",header=TRUE,sep="\t");
chrom_3R <- chrom_3R[which(chrom_3R$chromosome != 4),];
chrom_3R$chromosome <- factor(chrom_3R$chromosome,levels=c("X","2L","2R","3L","3R"));

chrom_X <- read.table("/Users/kraigrs/Wittkopp/gene_regulatory_network/Kellis_chromosomal_targets_X.txt",header=TRUE,sep="\t");
chrom_X <- chrom_X[which(chrom_X$chromosome != 4),];
chrom_X$chromosome <- factor(chrom_X$chromosome,levels=c("X","2L","2R","3L","3R"));

pdf(file="/Users/kraigrs/Desktop/Kellis_chromosomal_targets.pdf",height=8.5,width=11);

par(mfrow=c(2,3),oma = c(0,0,3,0));

boxplot(chr2L_regs/regulators ~ chromosome,data=chrom_2L,varwidth=TRUE,ylim=c(0,1),xlab="",ylab="proportion regulators",main="2L");
boxplot(chr3L_regs/regulators ~ chromosome,data=chrom_3L,varwidth=TRUE,ylim=c(0,1),xlab="",ylab="proportion regulators",main="3L");
boxplot(chrX_regs/regulators ~ chromosome,data=chrom_X,varwidth=TRUE,ylim=c(0,1),xlab="chromosome",ylab="proportion regulators",main="X");
boxplot(chr2R_regs/regulators ~ chromosome,data=chrom_2R,varwidth=TRUE,ylim=c(0,1),xlab="chromosome",ylab="proportion regulators",main="2R");
boxplot(chr3R_regs/regulators ~ chromosome,data=chrom_3R,varwidth=TRUE,ylim=c(0,1),xlab="chromosome",ylab="",main="3R");

mtext("Kellis data", outer = TRUE, cex = 1.5);

dev.off();