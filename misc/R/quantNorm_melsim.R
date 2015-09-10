source("/Users/kraigrs/Wittkopp/Rscripts/quantNorm.r");

mel_sim_par <- read.table("/Users/kraigrs/Wittkopp/mel_sim_data/mRNA-Seq/zhr+tsimbazaza/zhr+tsimbazaza.mosaik.exons.txt",header=T,sep="\t");
mel_sim_hyb1 <- read.table("/Users/kraigrs/Wittkopp/mel_sim_data/mRNA-Seq/zhrXtsimbazaza/zhrXtsimbazaza.mosaik.exons.txt",header=T,sep="\t");
mel_sim_hyb2 <- read.table("/Users/kraigrs/Wittkopp/mel_sim_data/mRNA-Seq/tsimbazazaXzhr/tsimbazazaXzhr.mosaik.exons.txt",header=T,sep="\t");

constExons <- read.table("/Users/kraigrs/Wittkopp/McManus/constitutive_regions.bed",sep="\t",skip=1);
lengths <- constExons[,3]-constExons[,2];
names <- paste(constExons[,4],"_",constExons[,2],",",constExons[,3],sep="");
constExons <- cbind(constExons,names,lengths);

mel_sim_parHyb1 <- merge(constExons,merge(mel_sim_par,mel_sim_hyb1,by.x="gene_exon",by.y="gene_exon"),by.x="names",by.y="gene_exon");
m <- cbind(mel_sim_parHyb1[,9],mel_sim_parHyb1[,10],mel_sim_parHyb1[,14],mel_sim_parHyb1[,15]);
mel_sim_parHyb1_qnorm <- quantNorm(m);

#see distribution of read count change
changes <- NULL;
for(i in 1:nrow(m))
{
  for(j in 1:ncol(m))
  {
    changes <- c(changes,m[i,j]-mel_sim_parHyb1_qnorm[i,j]);
  }
}
hist(changes);

changes <- NULL;
orig <- NULL;
for(i in 1:nrow(m))
{
  for(j in 1:ncol(m))
  {
    changes <- c(changes,abs(m[i,j]-mel_sim_parHyb1_qnorm[i,j]));
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
legend("topright",legend=c("zhr+tsim_zhr","zhr+tsim_tsim","zhrXtsim_zhr","zhrXtsim_tsim"),fill=c("blue","red","orange","green"));

# cdfs
plot(ecdf(log10(m[,1])),col="blue",xlim=c(0,6),ylim=c(0,1),xlab="",ylab="",main="",xaxt="n",yaxt="n");
par(new=TRUE);
plot(ecdf(log10(m[,2])),col="red",xlim=c(0,6),ylim=c(0,1),xlab="",ylab="",main="",xaxt="n",yaxt="n");
par(new=TRUE);
plot(ecdf(log10(m[,3])),col="orange",xlim=c(0,6),ylim=c(0,1),xlab="",ylab="",main="",xaxt="n",yaxt="n");
par(new=TRUE);
plot(ecdf(log10(m[,4])),col="green",xlim=c(0,6),ylim=c(0,1),xlab="log10(counts)",ylab="density",main="Distribution of counts per exon (before quantile normalization)");
legend("topright",legend=c("zhr+tsim_zhr","zhr+tsim_tsim","zhrXtsim_zhr","zhrXtsim_tsim"),fill=c("blue","red","orange","green"));

#compare log2 ratios before and after normalization
par(mfrow=c(1,2));
plot(log2(m[,1]/m[,2]),log2(floor(mel_sim_parHyb1_qnorm[,1])/floor(mel_sim_parHyb1_qnorm[,2])),col=rgb(0,0,0,0.2),cex=0.2,pch=19,xlab="log2(mel/sim) raw",ylab="log2(mel/sim) quantile normalized",main="zhr+tsimbazaza");
abline(a=0,b=1);
plot(log2(m[,3]/m[,4]),log2(floor(mel_sim_parHyb1_qnorm[,3])/floor(mel_sim_parHyb1_qnorm[,4])),col=rgb(0,0,0,0.2),cex=0.2,pch=19,xlab="log2(mel/sim) raw",ylab="log2(mel/sim) quantile normalized",main="zhrXtsimbazaza");
abline(a=0,b=1);

#floored
plot(log2(m[,1]/m[,2]),log2(floor(mel_sim_parHyb1_qnorm[,1])/floor(mel_sim_parHyb1_qnorm[,2])),col=rgb(0,0,0,0.2),cex=0.2,pch=19,xlab="log2(mel/sim) raw",ylab="log2(mel/sim) quantile normalized",main="zhr+tsimbazaza");
abline(a=0,b=1);

#floored
plot(log2(m[,3]/m[,4]),log2(floor(mel_sim_parHyb1_qnorm[,3])/floor(mel_sim_parHyb1_qnorm[,4])),col=rgb(0,0,0,0.2),cex=0.2,pch=19,xlab="log2(mel/sim) raw",ylab="log2(mel/sim) quantile normalized",main="zhrXtsimbazaza");
abline(a=0,b=1);