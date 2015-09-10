############################################################################################
# script to explore different properties of the data
############################################################################################

# Dpse_TL gDNA

Dpse_TL <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/Dpse_TL.genes.txt",header=TRUE,sep="\t");
prop <- round(sum(Dpse_TL$Dpse_TL)/(sum(Dpse_TL$Dpse_TL)+sum(Dpse_TL$Dbog_Toro1)),digits=3);

pdf("/Users/kraigrs/Wittkopp/Machado_Dpse/plots/Dpse_TL_gDNA.genes.pdf");
hist(Dpse_TL$Dpse_TL/(Dpse_TL$Dpse_TL+Dpse_TL$Dbog_Toro1),breaks=20,xlab="proportion of Dpse_TL allele",main="Dpse_TL_gDNA.genes");
abline(v=prop,col="red",lty=2);
legend("topleft",legend=paste("global proportion = ",prop,sep=""),lty=2,col="red");
dev.off();

# Dbog_Toro1 gDNA

Dbog_Toro1 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/Dbog_Toro1.genes.txt",header=TRUE,sep="\t");
prop <- round(sum(Dbog_Toro1$Dbog_Toro1)/(sum(Dbog_Toro1$Dpse_TL)+sum(Dbog_Toro1$Dbog_Toro1)),digits=3);

pdf("/Users/kraigrs/Wittkopp/Machado_Dpse/plots/Dbog_Toro1_gDNA.genes.pdf");
hist(Dbog_Toro1$Dbog_Toro1/(Dbog_Toro1$Dpse_TL+Dbog_Toro1$Dbog_Toro1),breaks=20,xlab="proportion of Dbog_Toro1 allele",main="Dbog_Toro1_gDNA.genes");
abline(v=prop,col="red",lty=2);
legend("topleft",legend=paste("global proportion = ",prop,sep=""),lty=2,col="red");
dev.off();

###################

Dpse_TL <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/Dpse_TL.genes.txt",header=TRUE,sep="\t");
Dbog_Toro1 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/Dbog_Toro1.genes.txt",header=TRUE,sep="\t");

data <- merge(Dpse_TL,Dbog_Toro1,by.x="gene",by.y="gene");

annotation <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/genomes/dpse_r3.1.genes.bed",header=FALSE,sep="\t");
colnames(annotation) <- c("chromosome","start","stop","gene","empty","strand");

seq_div <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/genomes/dpse_r3.1.Dpse_TL.Dbog_Toro1.BQSR_genes.seq_div.txt",header=TRUE,sep="\t");

info <- merge(annotation,seq_div,by.x="gene",by.y="gene");

data_info <- merge(data,info,by.x="gene",by.y="gene");
data_annotation <- subset(data_info,Dpse_TL.x+Dbog_Toro1.y>=20);
data_no_AS <- subset(data_info,Dpse_TL.x+Dbog_Toro1.y<20);

###############################
# make plots by chromosome

par(mfrow=c(1,2));

plot(1-data_annotation$Dpse_TL.x/(data_annotation$Dpse_TL.x+data_annotation$Dbog_Toro1.x),100*data_annotation$SNPs/data_annotation$coding,pch=19,cex=0.5,col=rgb(0,0,0,0.3),xlab="Dpse mismapping",ylab="% sequence divergence");

plot(1-data_annotation$Dbog_Toro1.y/(data_annotation$Dpse_TL.y+data_annotation$Dbog_Toro1.y),100*data_annotation$SNPs/data_annotation$coding,pch=19,cex=0.5,col=rgb(0,0,0,0.3),xlab="Dbog mismapping",ylab="");

data_XL <- data_annotation[grep("^XL",data_annotation$chromosome,perl=TRUE),];
data_4 <- data_annotation[grep("^4",data_annotation$chromosome,perl=TRUE),];
data_3 <- data_annotation[grep("^3",data_annotation$chromosome,perl=TRUE),];
data_XR <- data_annotation[grep("^XR",data_annotation$chromosome,perl=TRUE),];
data_2 <- data_annotation[grep("^2",data_annotation$chromosome,perl=TRUE),];

boxplot((data_XL$SNPs/data_XL$coding)*100,at=0.5,xlim=c(0,3),ylim=c(0,2),xlab="",ylab="% sequence divergence",main="Dpse and Dbog",outline=FALSE);
boxplot((data_4$SNPs/data_4$coding)*100,add=TRUE,at=1,xlim=c(0,3),ylim=c(0,2),xlab="",ylab="",main="",yaxt="n",outline=FALSE);
boxplot((data_3$SNPs/data_3$coding)*100,add=TRUE,at=1.5,xlim=c(0,3),ylim=c(0,2),xlab="",ylab="",main="",yaxt="n",outline=FALSE);
boxplot((data_XR$SNPs/data_XR$coding)*100,add=TRUE,at=2,xlim=c(0,3),ylim=c(0,2),xlab="",ylab="",main="",yaxt="n",outline=FALSE);
boxplot((data_2$SNPs/data_2$coding)*100,add=TRUE,at=2.5,xlim=c(0,3),ylim=c(0,2),xlab="",ylab="",main="",yaxt="n",outline=FALSE);
axis(side=1,at=c(0.5,1,1.5,2,2.5),labels=c("X","4","3","neo-X","2"));

par(mfrow=c(1,2));

boxplot(1-data_XL$Dpse_TL.x/(data_XL$Dpse_TL.x+data_XL$Dbog_Toro1.x),at=0.5,xlim=c(0,3),ylim=c(0,1),xlab="",ylab="proportion mismapping",main="Dpse",pch=16,cex=0.5);
boxplot(1-data_4$Dpse_TL.x/(data_4$Dpse_TL.x+data_4$Dbog_Toro1.x),add=TRUE,at=1,xlim=c(0,3),ylim=c(0,1),xlab="",ylab="",main="",yaxt="n",pch=16,cex=0.5);
boxplot(1-data_3$Dpse_TL.x/(data_3$Dpse_TL.x+data_3$Dbog_Toro1.x),add=TRUE,at=1.5,xlim=c(0,3),ylim=c(0,1),xlab="",ylab="",main="",yaxt="n",pch=16,cex=0.5);
boxplot(1-data_XR$Dpse_TL.x/(data_XR$Dpse_TL.x+data_XR$Dbog_Toro1.x),add=TRUE,at=2,xlim=c(0,3),ylim=c(0,1),xlab="",ylab="",main="",yaxt="n",pch=16,cex=0.5);
boxplot(1-data_2$Dpse_TL.x/(data_2$Dpse_TL.x+data_2$Dbog_Toro1.x),add=TRUE,at=2.5,xlim=c(0,3),ylim=c(0,1),xlab="",ylab="",main="",yaxt="n",pch=16,cex=0.5);
axis(side=1,at=c(0.5,1,1.5,2,2.5),labels=c("X","4","3","neo-X","2"));

boxplot(1-data_XL$Dbog_Toro1.y/(data_XL$Dpse_TL.y+data_XL$Dbog_Toro1.y),at=0.5,xlim=c(0,3),ylim=c(0,1),xlab="",ylab="",main="Dbog",pch=16,cex=0.5);
boxplot(1-data_4$Dbog_Toro1.y/(data_4$Dpse_TL.y+data_4$Dbog_Toro1.y),add=TRUE,at=1,xlim=c(0,3),ylim=c(0,1),xlab="",ylab="",main="",yaxt="n",pch=16,cex=0.5);
boxplot(1-data_3$Dbog_Toro1.y/(data_3$Dpse_TL.y+data_3$Dbog_Toro1.y),add=TRUE,at=1.5,xlim=c(0,3),ylim=c(0,1),xlab="",ylab="",main="",yaxt="n",pch=16,cex=0.5);
boxplot(1-data_XR$Dbog_Toro1.y/(data_XR$Dpse_TL.y+data_XR$Dbog_Toro1.y),add=TRUE,at=2,xlim=c(0,3),ylim=c(0,1),xlab="",ylab="",main="",yaxt="n",pch=16,cex=0.5);
boxplot(1-data_2$Dbog_Toro1.y/(data_2$Dpse_TL.y+data_2$Dbog_Toro1.y),add=TRUE,at=2.5,xlim=c(0,3),ylim=c(0,1),xlab="",ylab="",main="",yaxt="n",pch=16,cex=0.5);
axis(side=1,at=c(0.5,1,1.5,2,2.5),labels=c("X","4","3","neo-X","2"));

###############################

par(mfrow=c(1,2));
hist(log2(data_annotation$SNPs+1),breaks=100,xlab="log2(#SNPs+1)",main=">= 20 AS reads in Dpse and Dbog");
hist(log2(data_no_AS$SNPs+1),breaks=100,xlab="log2(#SNPs+1)",main="< 20 AS reads in Dpse or Dbog");

prop1 <- sum(data_annotation$Dpse_TL.x)/(sum(data_annotation$Dpse_TL.x)+sum(data_annotation$Dbog_Toro1.x));
prop2 <- sum(data_annotation$Dbog_Toro1.y)/(sum(data_annotation$Dpse_TL.y)+sum(data_annotation$Dbog_Toro1.y));

par(mfrow=c(1,2));
hist(Dpse_TL$Dpse_TL/(Dpse_TL$Dpse_TL+Dpse_TL$Dbog_Toro1),breaks=20,xlab="proportion of Dpse_TL allele",main="Dpse_TL_gDNA.genes");
abline(v=prop1,col="red",lty=2);
legend("topleft",legend=paste("global proportion = ",round(prop1,digits=3),sep=""),lty=2,col="red");
hist(Dbog_Toro1$Dbog_Toro1/(Dbog_Toro1$Dpse_TL+Dbog_Toro1$Dbog_Toro1),breaks=20,xlab="proportion of Dbog_Toro1 allele",main="Dbog_Toro1_gDNA.genes");
abline(v=prop2,col="red",lty=2);
legend("topleft",legend=paste("global proportion = ",round(prop2,digits=3),sep=""),lty=2,col="red");

nrow(Dpse_TL);
nrow(Dbog_Toro1);
nrow(data);
nrow(data_annotation);
nrow(subset(data,Dpse_TL.x/(Dpse_TL.x+Dbog_Toro1.x) >= 0.9 & Dbog_Toro1.y/(Dpse_TL.y+Dbog_Toro1.y) >= 0.9));
nrow(subset(data,Dpse_TL.x/(Dpse_TL.x+Dbog_Toro1.x) == 1 & Dbog_Toro1.y/(Dpse_TL.y+Dbog_Toro1.y) == 1));

temp <- subset(data_annotation,Dpse_TL.x/(Dpse_TL.x+Dbog_Toro1.x) >= 0.99 & Dbog_Toro1.y/(Dpse_TL.y+Dbog_Toro1.y) >= 0.99);
temp1 <- data.frame(temp[,1]);
colnames(temp1) <- "gene";

write.table(temp1,file="/Users/kraigrs/Wittkopp/Machado_Dpse/data/reliable_genes.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE);

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

plot(data_XL$genomic/1000000,1-data_XL$Dpse_TL.x/(data_XL$Dpse_TL.x+data_XL$Dbog_Toro1.x),xlab="",ylab="proportion Dpse_TL mismapping",pch=19,col="black",ylim=c(0,1),cex=0.5,main="",xlim=c(0.02,125),xaxt="n");
x <- data_XL$genomic/1000000;
y <- 1-data_XL$Dpse_TL.x/(data_XL$Dpse_TL.x+data_XL$Dbog_Toro1.x);
temp <- data.frame(x,y);
lo <- loess(y~x,data=temp);
lines(x[order(x)],lo$fitted[order(x)], col="red",type="l",lwd=3);

points(data_4$genomic/1000000,1-data_4$Dpse_TL.x/(data_4$Dpse_TL.x+data_4$Dbog_Toro1.x),col="gray",cex=0.5,pch=19);
x <- data_4$genomic/1000000;
y <- 1-data_4$Dpse_TL.x/(data_4$Dpse_TL.x+data_4$Dbog_Toro1.x);
temp <- data.frame(x,y);
lo <- loess(y~x,data=temp);
lines(x[order(x)],lo$fitted[order(x)], col="red",type="l",lwd=3);

points(data_3$genomic/1000000,1-data_3$Dpse_TL.x/(data_3$Dpse_TL.x+data_3$Dbog_Toro1.x),col="black",cex=0.5,pch=19);
x <- data_3$genomic/1000000;
y <- 1-data_3$Dpse_TL.x/(data_3$Dpse_TL.x+data_3$Dbog_Toro1.x);
temp <- data.frame(x,y);
lo <- loess(y~x,data=temp);
lines(x[order(x)],lo$fitted[order(x)], col="red",type="l",lwd=3);

points(data_XR$genomic/1000000,1-data_XR$Dpse_TL.x/(data_XR$Dpse_TL.x+data_XR$Dbog_Toro1.x),col="gray",cex=0.5,pch=19);
x <- data_XR$genomic/1000000;
y <- 1-data_XR$Dpse_TL.x/(data_XR$Dpse_TL.x+data_XR$Dbog_Toro1.x);
temp <- data.frame(x,y);
lo <- loess(y~x,data=temp);
lines(x[order(x)],lo$fitted[order(x)], col="red",type="l",lwd=3);

points(data_2$genomic/1000000,1-data_2$Dpse_TL.x/(data_2$Dpse_TL.x+data_2$Dbog_Toro1.x),col="black",cex=0.5,pch=19);
x <- data_2$genomic/1000000;
y <- 1-data_2$Dpse_TL.x/(data_2$Dpse_TL.x+data_2$Dbog_Toro1.x);
temp <- data.frame(x,y);
lo <- loess(y~x,data=temp);
lines(x[order(x)],lo$fitted[order(x)], col="red",type="l",lwd=3);



plot(data_XL$genomic/1000000,data_XL$Dpse_TL.y/(data_XL$Dpse_TL.y+data_XL$Dbog_Toro1.y),xlab="",ylab="proportion Dbog_Toro1 mismapping",pch=19,col="black",ylim=c(0,1),cex=0.5,main="",xlim=c(0.02,125),xaxt="n");
x <- data_XL$genomic/1000000;
y <- data_XL$Dpse_TL.y/(data_XL$Dpse_TL.y+data_XL$Dbog_Toro1.y);
temp <- data.frame(x,y);
lo <- loess(y~x,data=temp);
lines(x[order(x)],lo$fitted[order(x)], col="red",type="l",lwd=3);

points(data_4$genomic/1000000,data_4$Dpse_TL.y/(data_4$Dpse_TL.y+data_4$Dbog_Toro1.y),col="gray",cex=0.5,pch=19);
x <- data_4$genomic/1000000;
y <- data_4$Dpse_TL.y/(data_4$Dpse_TL.y+data_4$Dbog_Toro1.y);
temp <- data.frame(x,y);
lo <- loess(y~x,data=temp);
lines(x[order(x)],lo$fitted[order(x)], col="red",type="l",lwd=3);

points(data_3$genomic/1000000,data_3$Dpse_TL.y/(data_3$Dpse_TL.y+data_3$Dbog_Toro1.y),col="black",cex=0.5,pch=19);
x <- data_3$genomic/1000000;
y <- data_3$Dpse_TL.y/(data_3$Dpse_TL.y+data_3$Dbog_Toro1.y);
temp <- data.frame(x,y);
lo <- loess(y~x,data=temp);
lines(x[order(x)],lo$fitted[order(x)], col="red",type="l",lwd=3);

points(data_XR$genomic/1000000,data_XR$Dpse_TL.y/(data_XR$Dpse_TL.y+data_XR$Dbog_Toro1.y),col="gray",cex=0.5,pch=19);
x <- data_XR$genomic/1000000;
y <- data_XR$Dpse_TL.y/(data_XR$Dpse_TL.y+data_XR$Dbog_Toro1.y);
temp <- data.frame(x,y);
lo <- loess(y~x,data=temp);
lines(x[order(x)],lo$fitted[order(x)], col="red",type="l",lwd=3);

points(data_2$genomic/1000000,data_2$Dpse_TL.y/(data_2$Dpse_TL.y+data_2$Dbog_Toro1.y),col="black",cex=0.5,pch=19);
x <- data_2$genomic/1000000;
y <- data_2$Dpse_TL.y/(data_2$Dpse_TL.y+data_2$Dbog_Toro1.y);
temp <- data.frame(x,y);
lo <- loess(y~x,data=temp);
lines(x[order(x)],lo$fitted[order(x)], col="red",type="l",lwd=3);


par(mai = c(0.6, 0.6, 0.2, 0.2));

#plot(data_XL$genomic/1000000,100*(data_XL$SNPs/data_XL$coding),xlab="genomic position (Mb)",ylab="% sequence divergence",pch=19,col="black",ylim=c(0,2.5),cex=0.5,main="",xlim=c(0.02,125));
plot(data_XL$genomic/1000000,100*(data_XL$SNPs/data_XL$coding),xlab="",xaxt="n",ylab="% sequence divergence",pch=19,col="black",ylim=c(0,6.5),cex=0.5,main="",xlim=c(0.02,125));

locations <- c(((max(data_XL$genomic)+min(data_XL$genomic))/2)/1000000,
((max(data_4$genomic)+min(data_4$genomic))/2)/1000000,
((max(data_3$genomic)+min(data_3$genomic))/2)/1000000,
((max(data_XR$genomic)+min(data_XR$genomic))/2)/1000000,
((max(data_2$genomic)+min(data_2$genomic))/2)/1000000);

axis(side=1,at=locations,labels=c("X","4","3","neo-X","2"));

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




boxplot((data_XL$SNPs/data_XL$coding)*100,at=0.5,xlim=c(0,3),ylim=c(0,6.5),xlab="",ylab="",main="",yaxt="n",pch=16);
boxplot((data_4$SNPs/data_4$coding)*100,add=TRUE,at=1,xlim=c(0,3),ylim=c(0,6.5),xlab="",ylab="",main="",yaxt="n",pch=16);
boxplot((data_3$SNPs/data_3$coding)*100,add=TRUE,at=1.5,xlim=c(0,3),ylim=c(0,6.5),xlab="",ylab="",main="",yaxt="n",pch=16);
boxplot((data_XR$SNPs/data_XR$coding)*100,add=TRUE,at=2,xlim=c(0,3),ylim=c(0,6.5),xlab="",ylab="",main="",yaxt="n",pch=16);
boxplot((data_2$SNPs/data_2$coding)*100,add=TRUE,at=2.5,xlim=c(0,3),ylim=c(0,6.5),xlab="",ylab="",main="",yaxt="n",pch=16);
axis(side=1,at=c(0.5,1,1.5,2,2.5),labels=c("X","4","3","neo-X","2"));



