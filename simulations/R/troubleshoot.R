SNPs <- read.table("/Users/kraigrs/Wittkopp/Simulations/chr2L_SNPs_in_windows.txt",header=FALSE,sep="\t");
reads <- read.table("/Users/kraigrs/Wittkopp/Simulations/chr2L_reads_in_windows.txt",header=FALSE,sep="\t");

plot(log10(reads[,4]),SNPs[,4],xlab="log10(reads)",ylab="SNPs",main="Windows on chr2L",pch=19,col=rgb(0,0,0,0.5),cex=0.5);

plot(reads[,4],SNPs[,4],xlab="reads",ylab="SNPs",main="Windows on chr2L",pch=19,col=rgb(0,0,0,0.5),cex=0.5);



SNPs <- read.table("/Users/kraigrs/Wittkopp/Simulations/chr2L_SNPs_in_regs.txt",header=FALSE,sep="\t");
reads <- read.table("/Users/kraigrs/Wittkopp/Simulations/chr2L_reads_in_regs.txt",header=FALSE,sep="\t");
data <- merge(SNPs,reads,by.x="V4",by.y="V4");

expected <- 4*(reads[,3]-reads[,2]+1);

plot(data[,13]/(data[,4]-data[,3]),data[,7]/(data[,4]-data[,3]),xlab="reads",ylab="SNPs",main="Exons on chr2L",pch=19,col=rgb(0,0,0,0.5),cex=0.5);
