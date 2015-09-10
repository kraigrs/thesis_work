############################################################################################
# script to compare sequence divergence and regulatory divergence
############################################################################################

FET <- function(locus)
{
	table <- matrix(c(locus[1],locus[3],locus[2],locus[4]),
					nrow=2,
					dimnames=list(c("H1","H2"),c("Allele1","Allele2")));
	test <- fisher.test(table,or=1,alternative="two.sided",conf.level=0.95);
	result <- c(test$estimate,test$conf.int[1],test$conf.int[2],test$p.value);
	return(result);
}

# define datasets

#reg_div1 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_F_carcass_reg_div.txt",header=TRUE,sep="\t");
#MOI1 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_F_carcass_MOI.txt",header=TRUE,sep="\t");
#data1 <- merge(reg_div1,MOI1,by.x=1,by.y=1);

#reg_div2 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H5_ovaries_reg_div.txt",header=TRUE,sep="\t");
#MOI2 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H5_ovaries_MOI.txt",header=TRUE,sep="\t");
#data2 <- merge(reg_div2,MOI2,by.x=1,by.y=1);

#data <- merge(data1,data2,by.x="gene",by.y="gene");

#imprinting_results <- t(apply(as.matrix(data[,c(4,5,29,30)]),1,FET));
#imprinting_qval <- p.adjust(imprinting_results[,4],method="fdr");

#data <- cbind(data,imprinting_qval);

reg_div <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_F_carcass_reg_div.txt",header=TRUE,sep="\t");
MOI <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_F_carcass_MOI.txt",header=TRUE,sep="\t");
data <- merge(reg_div,MOI,by.x=1,by.y=1);

annotation <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/genomes/dpse_r3.1.genes.bed",header=FALSE,sep="\t");
colnames(annotation) <- c("chromosome","start","stop","gene","empty","strand");

seq_div <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/genomes/dpse_r3.1.Dpse_TL.Dbog_Toro1.genes.seq_div.txt",header=TRUE,sep="\t");

info <- merge(annotation,seq_div,by.x="gene",by.y="gene");

data_annotation <- merge(data,info,by.x="gene",by.y="gene");

# females

# combine chromosomes from XL, 4, and XR

data_XL_group1a <- data_annotation[grep("^XL_group1a",data_annotation$chromosome,perl=TRUE),];
data_XL_group1e <- data_annotation[grep("^XL_group1e",data_annotation$chromosome,perl=TRUE),];
data_XL_group3a <- data_annotation[grep("^XL_group3a",data_annotation$chromosome,perl=TRUE),];
data_XL_group3b <- data_annotation[grep("^XL_group3b",data_annotation$chromosome,perl=TRUE),];
coord1 <- (data_XL_group1a$start+data_XL_group1a$stop)/2;
coord2 <- (data_XL_group1e$start+data_XL_group1e$stop)/2 + max(coord1);
coord3 <- (data_XL_group3a$start+data_XL_group3a$stop)/2 + max(coord2);
coord4 <- (data_XL_group3b$start+data_XL_group3b$stop)/2 + max(coord3);
coord <- c(coord1,coord2,coord3,coord4);
temp <- rbind(data_XL_group1a,data_XL_group1e,data_XL_group3a,data_XL_group3b);
data_XL <- data.frame(temp,coord);

data_4_group1 <- data_annotation[grep("^4_group1",data_annotation$chromosome,perl=TRUE),];
data_4_group2 <- data_annotation[grep("^4_group2",data_annotation$chromosome,perl=TRUE),];
data_4_group3 <- data_annotation[grep("^4_group3",data_annotation$chromosome,perl=TRUE),];
data_4_group4 <- data_annotation[grep("^4_group4",data_annotation$chromosome,perl=TRUE),];
coord1 <- (data_4_group1$start+data_4_group1$stop)/2;
coord2 <- (data_4_group2$start+data_4_group2$stop)/2 + max(coord1);
coord3 <- (data_4_group3$start+data_4_group3$stop)/2 + max(coord2);
coord4 <- (data_4_group4$start+data_4_group4$stop)/2 + max(coord3);
coord <- c(coord1,coord2,coord3,coord4);
temp <- rbind(data_4_group1,data_4_group2,data_4_group3,data_4_group4);
data_4 <- data.frame(temp,coord);

data_XR_group3a <- data_annotation[grep("^XR_group3a",data_annotation$chromosome,perl=TRUE),];
data_XR_group5 <- data_annotation[grep("^XR_group5",data_annotation$chromosome,perl=TRUE),];
data_XR_group6 <- data_annotation[grep("^XR_group6",data_annotation$chromosome,perl=TRUE),];
data_XR_group8 <- data_annotation[grep("^XR_group8",data_annotation$chromosome,perl=TRUE),];
coord1 <- (data_XR_group3a$start+data_XR_group3a$stop)/2;
coord2 <- (data_XR_group5$start+data_XR_group5$stop)/2 + max(coord1);
coord3 <- (data_XR_group6$start+data_XR_group6$stop)/2 + max(coord2);
coord4 <- (data_XR_group8$start+data_XR_group8$stop)/2 + max(coord3);
coord <- c(coord1,coord2,coord3,coord4);
temp <- rbind(data_XR_group3a,data_XR_group5,data_XR_group6,data_XR_group8);
data_XR <- data.frame(temp,coord);

data_2 <- data_annotation[grep("^2",data_annotation$chromosome,perl=TRUE),];
coord <- (data_2$start+data_2$stop)/2;
data_2 <- data.frame(data_2,coord);

data_3 <- data_annotation[grep("^3",data_annotation$chromosome,perl=TRUE),];
coord <- (data_3$start+data_3$stop)/2;
data_3 <- data.frame(data_3,coord);

#data_unknown <- data_annotation[grep("^Unknown",data_annotation$chromosome,perl=TRUE),];

genomic <- data_XL$coord;
data_XL <- data.frame(data_XL,genomic);

genomic <- data_4$coord + max(data_XL$genomic);
data_4 <- data.frame(data_4,genomic);

genomic <- data_3$coord + max(data_4$genomic);
data_3 <- data.frame(data_3,genomic);

genomic <- data_XR$coord + max(data_3$genomic);
data_XR <- data.frame(data_XR,genomic);

genomic <- data_2$coord + max(data_XR$genomic);
data_2 <- data.frame(data_2,genomic);

par(mfrow=c(3,1));
par(mai = c(0.03, 0.6, 0.2, 0.2));

plot(data_XL$genomic/1000000,abs(log2(data_XL$Dpse_TL.x/data_XL$Dbog_Toro1.x)),xlab="",ylab="|log2(P1/P2)|",pch=19,col="black",ylim=c(0,6),cex=0.5,main="ovaries",xlim=c(0.02,125),xaxt="n");
points(data_4$genomic/1000000,abs(log2(data_4$Dpse_TL.x/data_4$Dbog_Toro1.x)),col="gray",cex=0.5,pch=19);
points(data_3$genomic/1000000,abs(log2(data_3$Dpse_TL.x/data_3$Dbog_Toro1.x)),col="black",cex=0.5,pch=19);
points(data_XR$genomic/1000000,abs(log2(data_XR$Dpse_TL.x/data_XR$Dbog_Toro1.x)),col="gray",cex=0.5,pch=19);
points(data_2$genomic/1000000,abs(log2(data_2$Dpse_TL.x/data_2$Dbog_Toro1.x)),col="black",cex=0.5,pch=19);

plot(data_XL$genomic/1000000,abs(log2(data_XL$Dpse_TL.y/data_XL$Dbog_Toro1.y)),xlab="",ylab="|log2(A1/A2)| H6",pch=19,col="black",ylim=c(0,6),cex=0.5,main="",xlim=c(0.02,125),xaxt="n");
points(data_4$genomic/1000000,abs(log2(data_4$Dpse_TL.y/data_4$Dbog_Toro1.y)),col="gray",cex=0.5,pch=19);
points(data_3$genomic/1000000,abs(log2(data_3$Dpse_TL.y/data_3$Dbog_Toro1.y)),col="black",cex=0.5,pch=19);
points(data_XR$genomic/1000000,abs(log2(data_XR$Dpse_TL.y/data_XR$Dbog_Toro1.y)),col="gray",cex=0.5,pch=19);
points(data_2$genomic/1000000,abs(log2(data_2$Dpse_TL.y/data_2$Dbog_Toro1.y)),col="black",cex=0.5,pch=19);

par(mai = c(0.6, 0.6, 0.2, 0.2));

plot(data_XL$genomic/1000000,100*(data_XL$SNPs/data_XL$coding),xlab="genomic position (Mb)",ylab="% sequence divergence",pch=19,col="black",ylim=c(0,2.5),cex=0.5,main="",xlim=c(0.02,125));
x <- data_XL$genomic/1000000;
y <- 100*(data_XL$SNPs/data_XL$coding);
temp <- data.frame(x,y);
lo <- loess(y~x,data=temp);
lines(x[order(x)],lo$fitted[order(x)], col="red",type="l",lwd=3);

points(data_4$genomic/1000000,100*(data_4$SNPs/data_4$coding),col="gray",cex=0.5,pch=19);
x <- data_4$genomic/1000000;
y <- 100*(data_4$SNPs/data_4$coding);
temp <- data.frame(x,y);
lo <- loess(y~x,data=temp);
lines(x[order(x)],lo$fitted[order(x)], col="red",type="l",lwd=3);

points(data_3$genomic/1000000,100*(data_3$SNPs/data_3$coding),col="black",cex=0.5,pch=19);
x <- data_3$genomic/1000000;
y <- 100*(data_3$SNPs/data_3$coding);
temp <- data.frame(x,y);
lo <- loess(y~x,data=temp);
lines(x[order(x)],lo$fitted[order(x)], col="red",type="l",lwd=3);

points(data_XR$genomic/1000000,100*(data_XR$SNPs/data_XR$coding),col="gray",cex=0.5,pch=19);
x <- data_XR$genomic/1000000;
y <- 100*(data_XR$SNPs/data_XR$coding);
temp <- data.frame(x,y);
lo <- loess(y~x,data=temp);
lines(x[order(x)],lo$fitted[order(x)], col="red",type="l",lwd=3);

points(data_2$genomic/1000000,100*(data_2$SNPs/data_2$coding),col="black",cex=0.5,pch=19);
x <- data_2$genomic/1000000;
y <- 100*(data_2$SNPs/data_2$coding);
temp <- data.frame(x,y);
lo <- loess(y~x,data=temp);
lines(x[order(x)],lo$fitted[order(x)], col="red",type="l",lwd=3);

# males

data <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_testes_reg_div_no_X.txt",header=TRUE,sep="\t");

seq_div <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/genomes/dpse_r3.1.Dpse_TL.Dbog_Toro1.genes.seq_div.txt",header=TRUE,sep="\t");

data_annotation <- merge(data,seq_div,by.x="gene",by.y="gene");


# combine chromosome 4

data_4_group1 <- data_annotation[grep("^4_group1",data_annotation$chromosome,perl=TRUE),];
data_4_group2 <- data_annotation[grep("^4_group2",data_annotation$chromosome,perl=TRUE),];
data_4_group3 <- data_annotation[grep("^4_group3",data_annotation$chromosome,perl=TRUE),];
data_4_group4 <- data_annotation[grep("^4_group4",data_annotation$chromosome,perl=TRUE),];
coord1 <- (data_4_group1$start+data_4_group1$stop)/2;
coord2 <- (data_4_group2$start+data_4_group2$stop)/2 + max(coord1);
coord3 <- (data_4_group3$start+data_4_group3$stop)/2 + max(coord2);
coord4 <- (data_4_group4$start+data_4_group4$stop)/2 + max(coord3);
coord <- c(coord1,coord2,coord3,coord4);
temp <- rbind(data_4_group1,data_4_group2,data_4_group3,data_4_group4);
data_4 <- data.frame(temp,coord);

data_3 <- data_annotation[grep("^3",data_annotation$chromosome,perl=TRUE),];
coord <- (data_3$start+data_3$stop)/2;
data_3 <- data.frame(data_3,coord);

data_2 <- data_annotation[grep("^2",data_annotation$chromosome,perl=TRUE),];
coord <- (data_2$start+data_2$stop)/2;
data_2 <- data.frame(data_2,coord);

#data_unknown <- data_annotation[grep("^Unknown",data_annotation$chromosome,perl=TRUE),];

genomic <- data_4$coord + 24550000;
data_4 <- data.frame(data_4,genomic);

genomic <- data_3$coord + max(data_4$genomic);
data_3 <- data.frame(data_3,genomic);

genomic <- data_2$coord + max(data_3$genomic) + 24560000;
data_2 <- data.frame(data_2,genomic);

par(mfrow=c(3,1));
par(mai = c(0.03, 0.6, 0.2, 0.2));

plot(data_4$genomic/1000000,abs(log2(data_4$Dpse_TL.x/data_4$Dbog_Toro1.x)),xlab="",ylab="|log2(P1/P2)|",pch=19,col="gray",ylim=c(0,6),cex=0.5,main="testes",xlim=c(0.02,125),xaxt="n");
points(data_3$genomic/1000000,abs(log2(data_3$Dpse_TL.x/data_3$Dbog_Toro1.x)),col="black",cex=0.5,pch=19);
points(data_2$genomic/1000000,abs(log2(data_2$Dpse_TL.x/data_2$Dbog_Toro1.x)),col="black",cex=0.5,pch=19);

plot(data_4$genomic/1000000,abs(log2(data_4$Dpse_TL.y/data_4$Dbog_Toro1.y)),xlab="",ylab="|log2(A1/A2)| H6",pch=19,col="gray",ylim=c(0,6),cex=0.5,main="",xlim=c(0.02,125),xaxt="n");
points(data_3$genomic/1000000,abs(log2(data_3$Dpse_TL.y/data_3$Dbog_Toro1.y)),col="black",cex=0.5,pch=19);
points(data_2$genomic/1000000,abs(log2(data_2$Dpse_TL.y/data_2$Dbog_Toro1.y)),col="black",cex=0.5,pch=19);

par(mai = c(0.6, 0.6, 0.2, 0.2));

plot(data_4$genomic/1000000,100*(data_4$SNPs/data_4$coding),xlab="genomic position (Mb)",ylab="% sequence divergence",pch=19,col="gray",ylim=c(0,2.5),cex=0.5,main="",xlim=c(0.02,125));
x <- data_4$genomic/1000000;
y <- 100*(data_4$SNPs/data_4$coding);
temp <- data.frame(x,y);
lo <- loess(y~x,data=temp);
lines(x[order(x)],lo$fitted[order(x)], col="red",type="l",lwd=3);

points(data_3$genomic/1000000,100*(data_3$SNPs/data_3$coding),col="black",cex=0.5,pch=19);
x <- data_3$genomic/1000000;
y <- 100*(data_3$SNPs/data_3$coding);
temp <- data.frame(x,y);
lo <- loess(y~x,data=temp);
lines(x[order(x)],lo$fitted[order(x)], col="red",type="l",lwd=3);

points(data_2$genomic/1000000,100*(data_2$SNPs/data_2$coding),col="black",cex=0.5,pch=19);
x <- data_2$genomic/1000000;
y <- 100*(data_2$SNPs/data_2$coding);
temp <- data.frame(x,y);
lo <- loess(y~x,data=temp);
lines(x[order(x)],lo$fitted[order(x)], col="red",type="l",lwd=3);