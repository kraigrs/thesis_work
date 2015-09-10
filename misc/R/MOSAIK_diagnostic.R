data <- read.table("/Users/kraigrs/Wittkopp/Simulations/chr2L_counts_SNPs_expressed.bed",header=FALSE,sep="\t");

plot(data$V8,data$V7,xlab="#SNPs in exon",ylab="trancsript abundance for exon",pch=19,cex=0.4,col=rgb(0,0,0,0.2));
plot(data$V8,data$V7,xlab="#SNPs in exon",ylab="trancsript abundance for exon",pch=19,cex=0.4,col=rgb(0,0,0,0.2),xlim=c(0,50),ylim=c(0,10000));

cor(data$V8,data$V7)^2;

length <- data$V3 - data$V2;
SNPdensity <- data$V8/(length/1000);