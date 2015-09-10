############################################################################################
# script to calculate rho values
############################################################################################

boot <- function(mat)
{
	set1 <- sample(seq(1,nrow(mat),1),nrow(mat),replace=TRUE);
	set2 <- sample(seq(1,nrow(mat),1),nrow(mat),replace=TRUE);
	val1 <- 1-cor(log10(mat$Dpse_TL[set1]),log10(mat$Dbog_Toro1[set1]),use="pairwise.complete.obs",method="spearman");
	val2 <- 1-cor(log10(mat$Hyb_Dpse[set2]),log10(mat$Hyb_Dbog[set2]),use="pairwise.complete.obs",method="spearman");
	return(c(val1,val2));
}

#################

data <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_ovaries_reg_div.txt",header=TRUE,sep="\t");

annotation <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/genomes/dpse_r3.1.genes.bed",header=FALSE,sep="\t");
colnames(annotation) <- c("chromosome","start","stop","gene","empty","strand");

data_annotation <- merge(data,annotation,by.x="gene",by.y="gene");

data_XL <- data_annotation[grep("^XL",data_annotation$chromosome,perl=TRUE),];
data_4 <- data_annotation[grep("^4",data_annotation$chromosome,perl=TRUE),];
data_3 <- data_annotation[grep("^3",data_annotation$chromosome,perl=TRUE),];
data_XR <- data_annotation[grep("^XR",data_annotation$chromosome,perl=TRUE),];
data_2 <- data_annotation[grep("^2",data_annotation$chromosome,perl=TRUE),];
data_all <- data_annotation[grep("^[^X]",data_annotation$chromosome,perl=TRUE),];

# total parental and hybrid expression

cor(log10(data_XL$par1.x),log10(data_XL$par2.x),use="pairwise.complete.obs",method="spearman");
cor(log10(data_4$par1.x),log10(data_4$par2.x),use="pairwise.complete.obs",method="spearman");
cor(log10(data_3$par1.x),log10(data_3$par2.x),use="pairwise.complete.obs",method="spearman");
cor(log10(data_XR$par1.x),log10(data_XR$par2.x),use="pairwise.complete.obs",method="spearman");
cor(log10(data_2$par1.x),log10(data_2$par2.x),use="pairwise.complete.obs",method="spearman");

cor(log10(data_XL$hyb.x),log10(data_XL$hyb.y),use="pairwise.complete.obs",method="spearman");
cor(log10(data_4$hyb.x),log10(data_4$hyb.y),use="pairwise.complete.obs",method="spearman");
cor(log10(data_3$hyb.x),log10(data_3$hyb.y),use="pairwise.complete.obs",method="spearman");
cor(log10(data_XR$hyb.x),log10(data_XR$hyb.y),use="pairwise.complete.obs",method="spearman");
cor(log10(data_2$hyb.x),log10(data_2$hyb.y),use="pairwise.complete.obs",method="spearman");

chrXL <- 1-cor(log10(data_XL$par1.x),log10(data_XL$par2.x),use="pairwise.complete.obs",method="spearman");
chr4 <- 1-cor(log10(data_4$par1.x),log10(data_4$par2.x),use="pairwise.complete.obs",method="spearman");
chr3 <- 1-cor(log10(data_3$par1.x),log10(data_3$par2.x),use="pairwise.complete.obs",method="spearman");
chrXR <- 1-cor(log10(data_XR$par1.x),log10(data_XR$par2.x),use="pairwise.complete.obs",method="spearman");
chr2 <- 1-cor(log10(data_2$par1.x),log10(data_2$par2.x),use="pairwise.complete.obs",method="spearman");
overall <- 1-cor(log10(data_all$par1.x),log10(data_all$par2.x),use="pairwise.complete.obs",method="spearman");

plot(seq(0,6,1),rep(overall,7),type="l",lty=2,ylim=c(0,0.4),ylab="1-rho",xlab="",xaxt="n",main="parent1 vs. parent2");
axis(1,at=seq(1,6,1),labels=c("XL","4","3","XR","2"));
points(1,chrXL);
points(2,chr4);
points(3,chr3);
points(4,chrXR);
points(5,chr2);

plot(((data_XL$start+data_XL$stop)/2)/1000,log2(data_XL$par1.x/data_XL$par2.x),xlab="",ylab="log2(A1/A2)",pch=19,col=rgb(0,0,0,0.5),ylim=c(-8,8),main="XL");
plot(((data_4$start+data_4$stop)/2)/1000,log2(data_4$par1.x/data_4$par2.x),xlab="",ylab="log2(A1/A2)",pch=19,col=rgb(0,0,0,0.5),ylim=c(-8,8),main="4");
plot(((data_3$start+data_3$stop)/2)/1000,log2(data_3$par1.x/data_3$par2.x),xlab="",ylab="log2(A1/A2)",pch=19,col=rgb(0,0,0,0.5),ylim=c(-8,8),main="3");
plot(((data_XR$start+data_XR$stop)/2)/1000,log2(data_XR$par1.x/data_XR$par2.x),xlab="",ylab="log2(A1/A2)",pch=19,col=rgb(0,0,0,0.5),ylim=c(-8,8),main="XR");
plot(((data_2$start+data_2$stop)/2)/1000,log2(data_2$par1.x/data_2$par2.x),xlab="gene midpoint (kb)",ylab="log2(A1/A2)",pch=19,col=rgb(0,0,0,0.5),ylim=c(-8,8),main="2");

chrXL <- 1-cor(log2(data_XL$par1.x/data_XL$par2.x),log2(data_XL$Dpse_TL.x.x/data_XL$Dbog_Toro1.x.x),use="pairwise.complete.obs",method="spearman");
chr4 <- 1-cor(log10(data_4$par1.x),log10(data_4$par2.x),use="pairwise.complete.obs",method="spearman");
chr3 <- 1-cor(log10(data_3$par1.x),log10(data_3$par2.x),use="pairwise.complete.obs",method="spearman");
chrXR <- 1-cor(log10(data_XR$par1.x),log10(data_XR$par2.x),use="pairwise.complete.obs",method="spearman");
chr2 <- 1-cor(log10(data_2$par1.x),log10(data_2$par2.x),use="pairwise.complete.obs",method="spearman");
unknown <- 1-cor(log10(data_unknown$par1.x),log10(data_unknown$par2.x),use="pairwise.complete.obs",method="spearman");
overall <- 1-cor(log10(data_annotation$par1.x),log10(data_annotation$par2.x),use="pairwise.complete.obs",method="spearman");

# allele-specific parental and hybrid expression

data <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_F_carcass_reg_div.txt",header=TRUE,sep="\t");

list <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/genomes/dpse_r3.1.genes.bed",header=FALSE,sep="\t");
colnames(list) <- c("chromosome","start","stop","gene","empty","strand");

##### FBGs #####

#FBid2GLEANR <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/genomes/FBid2GLEANR.txt",header=TRUE,sep="\t");
#annotation <- merge(FBid2GLEANR,list,by.x="gene",by.y="gene");

#Zhang_FBG <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/Zhang_Dpse_FBG.txt",header=TRUE,sep="\t");

#temp <- merge(annotation,Zhang_FBG,by.x="GLEANR",by.y="GleanR.ID");

#data_annotation <- merge(data,temp,by.x="gene",by.y="gene");

#FBG <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/fold_change_gonads_FBG.txt",header=TRUE,sep="\t");
#temp <- merge(FBG,list,by.x="gene",by.y="gene");

#data_annotation <- merge(data,temp,by.x="gene",by.y="gene");

data_annotation <- merge(data,list,by.x="gene",by.y="gene");

data_XL <- data_annotation[grep("^XL",data_annotation$chromosome,perl=TRUE),];
data_4 <- data_annotation[grep("^4",data_annotation$chromosome,perl=TRUE),];
data_3 <- data_annotation[grep("^3",data_annotation$chromosome,perl=TRUE),];
data_XR <- data_annotation[grep("^XR",data_annotation$chromosome,perl=TRUE),];
data_2 <- data_annotation[grep("^2",data_annotation$chromosome,perl=TRUE),];
data_all <- data_annotation[grep("^[^X]",data_annotation$chromosome,perl=TRUE),];

# females

# H6
set.seed(12345); # F_carcass
#set.seed(67891); # ovaries

# H5
#set.seed(23456); # F_carcass
#set.seed(78912); # ovaries

# bootstrap

chr_XL_boot <- t(replicate(10000,boot(data_XL)));
chr_4_boot <- t(replicate(10000,boot(data_4)));
chr_3_boot <- t(replicate(10000,boot(data_3)));
chr_XR_boot <- t(replicate(10000,boot(data_XR)));
chr_2_boot <- t(replicate(10000,boot(data_2)));
chr_all_boot <- t(replicate(10000,boot(data_all)));

# parent1 to parent2

pdf(file="/Users/kraigrs/Wittkopp/Machado_Dpse/plots/TL_Toro1_H6_ovaries_tot_expr_div.pdf",height=5,width=5);

chrXL <- 1-cor(log10(data_XL$Dpse_TL),log10(data_XL$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");
chr4 <- 1-cor(log10(data_4$Dpse_TL),log10(data_4$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");
chr3 <- 1-cor(log10(data_3$Dpse_TL),log10(data_3$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");
chrXR <- 1-cor(log10(data_XR$Dpse_TL),log10(data_XR$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");
chr2 <- 1-cor(log10(data_2$Dpse_TL),log10(data_2$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");
overall <- 1-cor(log10(data_all$Dpse_TL),log10(data_all$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");

plot(seq(0,6,1),rep(overall,7),type="l",lty=1,ylim=c(0,0.2),ylab="1-rho",xlab="",xaxt="n",main="parent1 vs. parent2",col="red");
lines(c(0,6),rep(quantile(chr_all_boot[,1],probs=0.025),2),col="red",lty=2);
lines(c(0,6),rep(quantile(chr_all_boot[,1],probs=0.975),2),col="red",lty=2);

axis(1,at=seq(1,5,1),labels=c("X","4","3","neo-X","2"));

points(1,chrXL,pch=16);
points(c(1,1),quantile(chr_XL_boot[,1],probs=c(0.025,0.975)),pch="-");
lines(c(1,1),quantile(chr_XL_boot[,1],probs=c(0.025,0.975)),lty=1,col="black");

points(2,chr4,pch=16);
points(c(2,2),quantile(chr_4_boot[,1],probs=c(0.025,0.975)),pch="-");
lines(c(2,2),quantile(chr_4_boot[,1],probs=c(0.025,0.975)),lty=1,col="black");

points(3,chr3,pch=16);
points(c(3,3),quantile(chr_3_boot[,1],probs=c(0.025,0.975)),pch="-");
lines(c(3,3),quantile(chr_3_boot[,1],probs=c(0.025,0.975)),lty=1,col="black");

points(4,chrXR,pch=16);
points(c(4,4),quantile(chr_XR_boot[,1],probs=c(0.025,0.975)),pch="-");
lines(c(4,4),quantile(chr_XR_boot[,1],probs=c(0.025,0.975)),lty=1,col="black");

points(5,chr2,pch=16);
points(c(5,5),quantile(chr_2_boot[,1],probs=c(0.025,0.975)),pch="-");
lines(c(5,5),quantile(chr_2_boot[,1],probs=c(0.025,0.975)),lty=1,col="black");

dev.off();

# ASE in hybrid

pdf(file="/Users/kraigrs/Wittkopp/Machado_Dpse/plots/TL_Toro1_H6_ovaries_cis_expr_div.pdf",height=5,width=5);

chrXL <- 1-cor(log10(data_XL$Hyb_Dpse),log10(data_XL$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");
chr4 <- 1-cor(log10(data_4$Hyb_Dpse),log10(data_4$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");
chr3 <- 1-cor(log10(data_3$Hyb_Dpse),log10(data_3$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");
chrXR <- 1-cor(log10(data_XR$Hyb_Dpse),log10(data_XR$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");
chr2 <- 1-cor(log10(data_2$Hyb_Dpse),log10(data_2$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");
overall <- 1-cor(log10(data_all$Hyb_Dpse),log10(data_all$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");

plot(seq(0,6,1),rep(overall,7),type="l",lty=1,ylim=c(0,0.2),ylab="",xlab="",xaxt="n",main="hybrid allele1 vs. allele2",col="red");
lines(c(0,6),rep(quantile(chr_all_boot[,2],probs=0.025),2),col="red",lty=2);
lines(c(0,6),rep(quantile(chr_all_boot[,2],probs=0.975),2),col="red",lty=2);

axis(1,at=seq(1,5,1),labels=c("X","4","3","neo-X","2"));

points(1,chrXL,pch=16);
points(c(1,1),quantile(chr_XL_boot[,2],probs=c(0.025,0.975)),pch="-");
lines(c(1,1),quantile(chr_XL_boot[,2],probs=c(0.025,0.975)),lty=1,col="black");

points(2,chr4,pch=16);
points(c(2,2),quantile(chr_4_boot[,2],probs=c(0.025,0.975)),pch="-");
lines(c(2,2),quantile(chr_4_boot[,2],probs=c(0.025,0.975)),lty=1,col="black");

points(3,chr3,pch=16);
points(c(3,3),quantile(chr_3_boot[,2],probs=c(0.025,0.975)),pch="-");
lines(c(3,3),quantile(chr_3_boot[,2],probs=c(0.025,0.975)),lty=1,col="black");

points(4,chrXR,pch=16);
points(c(4,4),quantile(chr_XR_boot[,2],probs=c(0.025,0.975)),pch="-");
lines(c(4,4),quantile(chr_XR_boot[,2],probs=c(0.025,0.975)),lty=1,col="black");

points(5,chr2,pch=16);
points(c(5,5),quantile(chr_2_boot[,2],probs=c(0.025,0.975)),pch="-");
lines(c(5,5),quantile(chr_2_boot[,2],probs=c(0.025,0.975)),lty=1,col="black");

dev.off();

# males

data <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_M_carcass_reg_div.txt",header=TRUE,sep="\t");

list <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/genomes/dpse_r3.1.genes.bed",header=FALSE,sep="\t");
colnames(list) <- c("chromosome","start","stop","gene","empty","strand");

##### MBGs #####

#FBid2GLEANR <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/genomes/FBid2GLEANR.txt",header=TRUE,sep="\t");
#annotation <- merge(FBid2GLEANR,list,by.x="gene",by.y="gene");

#Zhang_MBG <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/Zhang_Dpse_MBG.txt",header=TRUE,sep="\t");

#temp <- merge(annotation,Zhang_MBG,by.x="GLEANR",by.y="GleanR.ID");

#data_annotation <- merge(data,temp,by.x="gene",by.y="gene");

#MBG <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/fold_change_gonads_MBG.txt",header=TRUE,sep="\t");
MBG <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/testis_specific.txt",header=TRUE,sep="\t");

temp <- merge(MBG,list,by.x="gene",by.y="gene");

data_annotation <- merge(data,temp,by.x="gene",by.y="gene");

################

#data_annotation <- merge(data,list,by.x="gene",by.y="gene");

data_XL <- data_annotation[grep("^XL",data_annotation$chromosome,perl=TRUE),];
data_4 <- data_annotation[grep("^4",data_annotation$chromosome,perl=TRUE),];
data_3 <- data_annotation[grep("^3",data_annotation$chromosome,perl=TRUE),];
data_XR <- data_annotation[grep("^XR",data_annotation$chromosome,perl=TRUE),];
data_2 <- data_annotation[grep("^2",data_annotation$chromosome,perl=TRUE),];
data_all <- data_annotation[grep("^[^X]",data_annotation$chromosome,perl=TRUE),];

#set.seed(34567); # H6 M_carcass
#set.seed(89123); # H6 testes
#set.seed(45678); # H5 M_carcass
#set.seed(91234); # H5 testes

# bootstrap

chr_XL_boot <- t(replicate(10000,boot(data_XL)));
chr_4_boot <- t(replicate(10000,boot(data_4)));
chr_3_boot <- t(replicate(10000,boot(data_3)));
chr_XR_boot <- t(replicate(10000,boot(data_XR)));
chr_2_boot <- t(replicate(10000,boot(data_2)));
chr_all_boot <- t(replicate(10000,boot(data_all)));

# parent1 to parent2

pdf(file="/Users/kraigrs/Wittkopp/Machado_Dpse/plots/TL_Toro1_H6_testes_tot_expr_div.pdf",height=5,width=5);

chrXL <- 1-cor(log10(data_XL$Dpse_TL),log10(data_XL$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");
chr4 <- 1-cor(log10(data_4$Dpse_TL),log10(data_4$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");
chr3 <- 1-cor(log10(data_3$Dpse_TL),log10(data_3$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");
chrXR <- 1-cor(log10(data_XR$Dpse_TL),log10(data_XR$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");
chr2 <- 1-cor(log10(data_2$Dpse_TL),log10(data_2$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");
overall <- 1-cor(log10(data_all$Dpse_TL),log10(data_all$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");

#par(mfrow=c(1,2));

#ylim=c(0,0.25),

plot(seq(0,5,1),rep(overall,6),type="l",lty=1,ylab="1-rho",xlab="",xaxt="n",main="parent1 vs. parent2",col="red");
lines(c(0,5),rep(quantile(chr_all_boot[,1],probs=0.025),2),col="red",lty=2);
lines(c(0,5),rep(quantile(chr_all_boot[,1],probs=0.975),2),col="red",lty=2);

axis(1,at=seq(0.5,5,1),labels=c("X","4","3","neo-X","2"));

points(0.5,chrXL,pch=16);
points(c(0.5,0.5),quantile(chr_XL_boot[,1],probs=c(0.025,0.975)),pch="-");
lines(c(0.5,0.5),quantile(chr_XL_boot[,1],probs=c(0.025,0.975)),lty=1,col="black");

points(1.5,chr4,pch=16);
points(c(1.5,1.5),quantile(chr_4_boot[,1],probs=c(0.025,0.975)),pch="-");
lines(c(1.5,1.5),quantile(chr_4_boot[,1],probs=c(0.025,0.975)),lty=1,col="black");

points(2.5,chr3,pch=16);
points(c(2.5,2.5),quantile(chr_3_boot[,1],probs=c(0.025,0.975)),pch="-");
lines(c(2.5,2.5),quantile(chr_3_boot[,1],probs=c(0.025,0.975)),lty=1,col="black");

points(3.5,chrXR,pch=16);
points(c(3.5,3.5),quantile(chr_XR_boot[,1],probs=c(0.025,0.975)),pch="-");
lines(c(3.5,3.5),quantile(chr_XR_boot[,1],probs=c(0.025,0.975)),lty=1,col="black");

points(4.5,chr2,pch=16);
points(c(4.5,4.5),quantile(chr_2_boot[,1],probs=c(0.025,0.975)),pch="-");
lines(c(4.5,4.5),quantile(chr_2_boot[,1],probs=c(0.025,0.975)),lty=1,col="black");

dev.off();

# ASE in H6

chr4 <- 1-cor(log10(data_4$Hyb_Dpse),log10(data_4$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");
chr3 <- 1-cor(log10(data_3$Hyb_Dpse),log10(data_3$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");
chr2 <- 1-cor(log10(data_2$Hyb_Dpse),log10(data_2$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");
overall <- 1-cor(log10(data_all$Hyb_Dpse),log10(data_all$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");

plot(seq(0,4,1),rep(overall,5),type="l",lty=1,ylim=c(0,0.4),ylab="1-rho",xlab="",xaxt="n",main="hybrid allele1 vs. allele2",col="red");
lines(c(0,4),rep(quantile(chr_all_boot[,2],probs=0.025),2),col="red",lty=2);
lines(c(0,4),rep(quantile(chr_all_boot[,2],probs=0.975),2),col="red",lty=2);

axis(1,at=seq(1,3,1),labels=c("4","3","2"));

points(1,chr4,pch=16);
points(c(1,1),quantile(chr_4_boot[,2],probs=c(0.025,0.975)),pch="-");
lines(c(1,1),quantile(chr_4_boot[,2],probs=c(0.025,0.975)),lty=1,col="black");

points(2,chr3,pch=16);
points(c(2,2),quantile(chr_3_boot[,2],probs=c(0.025,0.975)),pch="-");
lines(c(2,2),quantile(chr_3_boot[,2],probs=c(0.025,0.975)),lty=1,col="black");

points(3,chr2,pch=16);
points(c(3,3),quantile(chr_2_boot[,2],probs=c(0.025,0.975)),pch="-");
lines(c(3,3),quantile(chr_2_boot[,2],probs=c(0.025,0.975)),lty=1,col="black");






plot(((data_XL$start+data_XL$stop)/2)/1000,log2(data_XL$Dpse_TL.y.x/data_XL$Dbog_Toro1.y.x),xlab="",ylab="log2(A1/A2)",pch=19,col=rgb(0,0,0,0.5),ylim=c(-8,8));
plot(((data_4$start+data_4$stop)/2)/1000,log2(data_4$Dpse_TL.y.x/data_4$Dbog_Toro1.y.x),xlab="",ylab="log2(A1/A2)",pch=19,col=rgb(0,0,0,0.5),ylim=c(-8,8));
plot(((data_3$start+data_3$stop)/2)/1000,log2(data_3$Dpse_TL.y.x/data_3$Dbog_Toro1.y.x),xlab="",ylab="log2(A1/A2)",pch=19,col=rgb(0,0,0,0.5),ylim=c(-8,8));
plot(((data_XR$start+data_XR$stop)/2)/1000,log2(data_XR$Dpse_TL.y.x/data_XR$Dbog_Toro1.y.x),xlab="",ylab="log2(A1/A2)",pch=19,col=rgb(0,0,0,0.5),ylim=c(-8,8));
plot(((data_2$start+data_2$stop)/2)/1000,log2(data_2$Dpse_TL.y.x/data_2$Dbog_Toro1.y.x),xlab="gene midpoint (kb)",ylab="log2(A1/A2)",pch=19,col=rgb(0,0,0,0.5),ylim=c(-8,8));



par(mfrow=c(2,4));

rho <- round(cor(log10(data_annotation$par1.x),log10(data_annotation$hyb.x),use="pairwise.complete.obs",method="spearman"),digits=3);
plot(log10(data_annotation$par1.x),log10(data_annotation$hyb.x),xlim=c(1,6),ylim=c(1,6),xlab="log10(Dpse)",ylab="log10(Dpse X Dbog)",pch=19,cex=0.5,col=rgb(0,0,0,0.5));
legend("topright",legend=paste("rho = ",rho,sep=""));

rho <- round(cor(log10(data_annotation$par2.x),log10(data_annotation$hyb.x),use="pairwise.complete.obs",method="spearman"),digits=3);
plot(log10(data_annotation$par2.x),log10(data_annotation$hyb.x),xlim=c(1,6),ylim=c(1,6),xlab="log10(Dbog)",ylab="log10(Dpse X Dbog)",pch=19,cex=0.5,col=rgb(0,0,0,0.5));
legend("topright",legend=paste("rho = ",rho,sep=""));

rho <- round(cor(log10(data_annotation$par1.x),log10(data_annotation$hyb.y),use="pairwise.complete.obs",method="spearman"),digits=3);
plot(log10(data_annotation$par1.x),log10(data_annotation$hyb.y),xlim=c(1,6),ylim=c(1,6),xlab="log10(Dpse)",ylab="log10(Dbog X Dpse)",pch=19,cex=0.5,col=rgb(0,0,0,0.5));
legend("topright",legend=paste("rho = ",rho,sep=""));

rho <- round(cor(log10(data_annotation$par2.x),log10(data_annotation$hyb.y),use="pairwise.complete.obs",method="spearman"),digits=3);
plot(log10(data_annotation$par2.x),log10(data_annotation$hyb.y),xlim=c(1,6),ylim=c(1,6),xlab="log10(Dbog)",ylab="log10(Dbog X Dpse)",pch=19,cex=0.5,col=rgb(0,0,0,0.5));
legend("topright",legend=paste("rho = ",rho,sep=""));

rho <- round(cor(log10(data_annotation$par1.x),log10(data_annotation$par2.x),use="pairwise.complete.obs",method="spearman"),digits=3);
plot(log10(data_annotation$par1.x),log10(data_annotation$par2.x),xlim=c(1,6),ylim=c(1,6),xlab="log10(Dpse)",ylab="log10(Dbog)",pch=19,cex=0.5,col=rgb(0,0,0,0.5));
legend("topright",legend=paste("rho = ",rho,sep=""));

rho <- round(cor(log10(data_annotation$hyb.x),log10(data_annotation$hyb.y),use="pairwise.complete.obs",method="spearman"),digits=3);
plot(log10(data_annotation$hyb.x),log10(data_annotation$hyb.y),xlim=c(1,6),ylim=c(1,6),xlab="log10(Dpse X Dbog)",ylab="log10(Dbog X Dpse)",pch=19,cex=0.5,col=rgb(0,0,0,0.5));
legend("topright",legend=paste("rho = ",rho,sep=""));

plot(log2(data_annotation$Dpse_TL.y.x/data_annotation$Dbog_Toro1.y.x),log2(data_annotation$Dpse_TL.y.y/data_annotation$Dbog_Toro1.y.y),xlim=c(-10,10),ylim=c(-10,10),pch=19,cex=0.5,col=rgb(0,0,0,0.5),xlab="log2(Dpse/Dbog) H6",ylab="log2(Dpse/Dbog) H5");












