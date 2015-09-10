source("/Users/kraigrs/Wittkopp/Rscripts/quantNorm.r");

mel_mel_zhr <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/mRNA-Seq/zhr/zhr.mosaik.exons.txt",header=T,sep="\t");
mel_mel_z30 <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/mRNA-Seq/z30/z30.mosaik.exons.txt",header=T,sep="\t");
mel_mel_par <- merge(mel_mel_zhr[,c(1,2)],mel_mel_z30[,c(1,3)],by.x="gene_exon",by.y="gene_exon"); 

mel_mel_hyb1 <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mosaik.exons.txt",header=T,sep="\t");
mel_mel_hyb2 <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mosaik.exons.txt",header=T,sep="\t");

constExons <- read.table("/Users/kraigrs/Wittkopp/McManus/constitutive_regions.bed",sep="\t",skip=1);
lengths <- constExons[,3]-constExons[,2];
names <- paste(constExons[,4],"_",constExons[,2],",",constExons[,3],sep="");
constExons <- cbind(constExons,names,lengths);

mel_mel_parHyb1 <- merge(constExons,merge(mel_mel_par,mel_mel_hyb1,by.x="gene_exon",by.y="gene_exon"),by.x="names",by.y="gene_exon");
m1 <- cbind(mel_mel_parHyb1[,9],mel_mel_parHyb1[,10],mel_mel_parHyb1[,11],mel_mel_parHyb1[,12]);
mel_mel_parHyb1_qnorm <- quantNorm(m1);

mel_mel_parHyb2 <- merge(constExons,merge(mel_mel_par,mel_mel_hyb2,by.x="gene_exon",by.y="gene_exon"),by.x="names",by.y="gene_exon");
m2 <- cbind(mel_mel_parHyb2[,9],mel_mel_parHyb2[,10],mel_mel_parHyb2[,11],mel_mel_parHyb2[,12]);
mel_mel_parHyb2_qnorm <- quantNorm(m2);

mel_mel_hyb1Hyb2 <- merge(constExons,merge(mel_mel_hyb1,mel_mel_hyb2,by.x="gene_exon",by.y="gene_exon"),by.x="names",by.y="gene_exon");
m3 <- cbind(mel_mel_hyb1Hyb2[,9], mel_mel_hyb1Hyb2[,10], mel_mel_hyb1Hyb2[,14], mel_mel_hyb1Hyb2[,15]);
mel_mel_hyb1Hyb2_qnorm <- quantNorm(m3);

#cut 20 

mel_mel_parHyb1_cut20 <- subset(mel_mel_parHyb1,mel_mel_parHyb1[,9]+mel_mel_parHyb1[,10]>=20 & mel_mel_parHyb1[,11]+mel_mel_parHyb1[,12]>=20);
m1_cut20 <- cbind(mel_mel_parHyb1_cut20[,9],mel_mel_parHyb1_cut20[,10],mel_mel_parHyb1_cut20[,11],mel_mel_parHyb1_cut20[,12]);
mel_mel_parHyb1_cut20_qnorm <- quantNorm(m1_cut20);

mel_mel_parHyb2_cut20 <- subset(mel_mel_parHyb2,mel_mel_parHyb2[,9]+mel_mel_parHyb2[,10]>=20 & mel_mel_parHyb2[,11]+mel_mel_parHyb2[,12]>=20);
m2_cut20 <- cbind(mel_mel_parHyb2_cut20[,9],mel_mel_parHyb2_cut20[,10],mel_mel_parHyb2_cut20[,11],mel_mel_parHyb2_cut20[,12]);
mel_mel_parHyb2_cut20_qnorm <- quantNorm(m2_cut20);

mel_mel_hyb1Hyb2_cut20 <- subset(mel_mel_hyb1Hyb2, mel_mel_hyb1Hyb2[,9]+mel_mel_hyb1Hyb2[,10]>=20 & mel_mel_hyb1Hyb2[,14]+mel_mel_hyb1Hyb2[,15]>=20);
m3_cut20 <- cbind(mel_mel_hyb1Hyb2_cut20[,9],mel_mel_hyb1Hyb2_cut20[,10],mel_mel_hyb1Hyb2_cut20[,14],mel_mel_hyb1Hyb2_cut20[,15]);
mel_mel_hyb1Hyb2_cut20_qnorm <- quantNorm(m3_cut20);

#see distribution of read count change
changes <- NULL;
orig <- NULL;
for(i in 1:nrow(m))
{
  for(j in 1:ncol(m))
  {
    changes <- c(changes,m[i,j]-mel_mel_parHyb1_qnorm[i,j]);
    orig <- c(orig,m[i,j]);
  }
}
hist(changes);

changes <- NULL;
orig <- NULL;
for(i in 1:nrow(m))
{
  for(j in 1:ncol(m))
  {
    changes <- c(changes,abs(m[i,j]-mel_mel_parHyb1_qnorm[i,j]));
    orig <- c(orig,m[i,j]);
  }
}
plot(orig,changes,xlim=c(0,100000),ylim=c(-20000,20000),pch=19,cex=0.2,col=rgb(0,0,0,0.2));

# see distributions before and after
plot(density(log10(m[,1])),col="blue",xlim=c(0,6),ylim=c(0,0.4),xlab="",ylab="",main="",xaxt="n",yaxt="n");
par(new=TRUE);
plot(density(log10(m[,2])),col="red",xlim=c(0,6),ylim=c(0,0.4),xlab="",ylab="",main="",xaxt="n",yaxt="n");
par(new=TRUE);
plot(density(log10(m[,3])),col="orange",xlim=c(0,6),ylim=c(0,0.4),xlab="",ylab="",main="",xaxt="n",yaxt="n");
par(new=TRUE);
plot(density(log10(m[,4])),col="green",xlim=c(0,6),ylim=c(0,0.4),xlab="log10(counts)",ylab="density",main="Distribution of counts per exon (before quantile normalization)");
legend("topright",legend=c("zhr","z30","zhrXz30_zhr","zhrXz30_z30"),fill=c("blue","red","orange","green"));

# cdfs
plot(ecdf(log10(m[,1])),col="blue",xlim=c(0,6),ylim=c(0,1),xlab="",ylab="",main="",xaxt="n",yaxt="n");
par(new=TRUE);
plot(ecdf(log10(m[,2])),col="red",xlim=c(0,6),ylim=c(0,1),xlab="",ylab="",main="",xaxt="n",yaxt="n");
par(new=TRUE);
plot(ecdf(log10(m[,3])),col="orange",xlim=c(0,6),ylim=c(0,1),xlab="",ylab="",main="",xaxt="n",yaxt="n");
par(new=TRUE);
plot(ecdf(log10(m[,4])),col="green",xlim=c(0,6),ylim=c(0,1),xlab="log10(counts)",ylab="density",main="Distribution of counts per exon (before quantile normalization)");
legend("topright",legend=c("zhr","z30","zhrXz30_zhr","zhrXz30_z30"),fill=c("blue","red","orange","green"));

#compare log2 ratios before and after normalization
par(mfrow=c(1,2));
plot(log2(m[,1]/m[,2]),log2(floor(mel_mel_parHyb1_qnorm[,1])/floor(mel_mel_parHyb1_qnorm[,2])),col=rgb(0,0,0,0.2),cex=0.2,pch=19,xlab="log2(zhr/z30) raw",ylab="log2(zhr/z30) quantile normalized",main="zhr+z30");
abline(a=0,b=1);
plot(log2(m[,3]/m[,4]),log2(floor(mel_mel_parHyb1_qnorm[,3])/floor(mel_mel_parHyb1_qnorm[,4])),col=rgb(0,0,0,0.2),cex=0.2,pch=19,xlab="log2(zhr/z30) raw",ylab="log2(zhr/z30) quantile normalized",main="zhrXz30");
abline(a=0,b=1);

cols <- c("locus","ParS1","ParS2","ParBoth","Hyb1S1","Hyb1S2","Hyb1Both");

zhr_z30_parHyb1_qnorm <- cbind(names,floor(mel_mel_parHyb1_qnorm[,1]),floor(mel_mel_parHyb1_qnorm[,2]),rep(0,nrow(mel_mel_parHyb1_qnorm)),floor(mel_mel_parHyb1_qnorm[,3]),floor(mel_mel_parHyb1_qnorm[,4]),rep(0,nrow(mel_mel_parHyb1_qnorm)));
write.table(zhr_z30_parHyb1_qnorm,file="/Users/kraigrs/Wittkopp/mel_mel_data/zhr_z30_parHyb1_qnorm.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=cols);

cols <- c("locus","ParS1","ParS2","ParBoth","Hyb2S1","Hyb2S2","Hyb2Both");

zhr_z30_parHyb2_qnorm <- cbind(names,floor(mel_mel_parHyb2_qnorm[,1]),floor(mel_mel_parHyb2_qnorm[,2]),rep(0,nrow(mel_mel_parHyb2_qnorm)),floor(mel_mel_parHyb2_qnorm[,3]),floor(mel_mel_parHyb2_qnorm[,4]),rep(0,nrow(mel_mel_parHyb2_qnorm)));
write.table(zhr_z30_parHyb2_qnorm,file="/Users/kraigrs/Wittkopp/mel_mel_data/zhr_z30_parHyb2_qnorm.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=cols);

cols <- c("locus","Hyb1S1","Hyb1S2","Hyb1Both","Hyb2S1","Hyb2S2","Hyb2Both");

zhr_z30_hyb1Hyb2_qnorm <- cbind(names,floor(mel_mel_hyb1Hyb2_qnorm[,1]),floor(mel_mel_hyb1Hyb2_qnorm[,2]),rep(0,nrow(mel_mel_hyb1Hyb2_qnorm)),floor(mel_mel_hyb1Hyb2_qnorm[,3]),floor(mel_mel_hyb1Hyb2_qnorm[,4]),rep(0,nrow(mel_mel_hyb1Hyb2_qnorm)));
write.table(zhr_z30_hyb1Hyb2_qnorm,file="/Users/kraigrs/Wittkopp/mel_mel_data/zhr_z30_hyb1Hyb2_qnorm.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=cols);

cols <- c("locus","ParS1","ParS2","ParBoth","Hyb1S1","Hyb1S2","Hyb1Both");

names1 <- data.frame(mel_mel_parHyb1_cut20[,1]);
zhr_z30_parHyb1_cut20_qnorm <- cbind(names1,floor(mel_mel_parHyb1_cut20_qnorm[,1]),floor(mel_mel_parHyb1_cut20_qnorm[,2]),rep(0,nrow(mel_mel_parHyb1_cut20_qnorm)),floor(mel_mel_parHyb1_cut20_qnorm[,3]),floor(mel_mel_parHyb1_cut20_qnorm[,4]),rep(0,nrow(mel_mel_parHyb1_cut20_qnorm)));
write.table(zhr_z30_parHyb1_cut20_qnorm,file="/Users/kraigrs/Wittkopp/mel_mel_data/zhr_z30_parHyb1_cut20_qnorm.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=cols);

cols <- c("locus","ParS1","ParS2","ParBoth","Hyb2S1","Hyb2S2","Hyb2Both");

names2 <- data.frame(mel_mel_parHyb2_cut20[,1]);
zhr_z30_parHyb2_cut20_qnorm <- cbind(names2,floor(mel_mel_parHyb2_cut20_qnorm[,1]),floor(mel_mel_parHyb2_cut20_qnorm[,2]),rep(0,nrow(mel_mel_parHyb2_cut20_qnorm)),floor(mel_mel_parHyb2_cut20_qnorm[,3]),floor(mel_mel_parHyb2_cut20_qnorm[,4]),rep(0,nrow(mel_mel_parHyb2_cut20_qnorm)));
write.table(zhr_z30_parHyb2_cut20_qnorm,file="/Users/kraigrs/Wittkopp/mel_mel_data/zhr_z30_parHyb2_cut20_qnorm.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=cols);

cols <- c("locus","Hyb1S1","Hyb1S2","Hyb1Both","Hyb2S1","Hyb2S2","Hyb2Both");

names3 <- data.frame(mel_mel_hyb1Hyb2_cut20[,1]);
zhr_z30_hyb1Hyb2_cut20_qnorm <- cbind(names3,floor(mel_mel_hyb1Hyb2_cut20_qnorm[,1]),floor(mel_mel_hyb1Hyb2_cut20_qnorm[,2]),rep(0,nrow(mel_mel_hyb1Hyb2_cut20_qnorm)),floor(mel_mel_hyb1Hyb2_cut20_qnorm[,3]),floor(mel_mel_hyb1Hyb2_cut20_qnorm[,4]),rep(0,nrow(mel_mel_hyb1Hyb2_cut20_qnorm)));
write.table(zhr_z30_hyb1Hyb2_cut20_qnorm,file="/Users/kraigrs/Wittkopp/mel_mel_data/zhr_z30_hyb1Hyb2_cut20_qnorm.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=cols);










