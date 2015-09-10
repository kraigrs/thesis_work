path2data <- "/Users/wittkopp-lab/Desktop/Kraig/Wittkopp/Bowtie/mel_sim_data/Hybrids/s_2_sequence"

mate1 <- read.table(paste(path2data,"/s_2_sequence.mate1.avgQualPerBP.txt",sep=""),header=TRUE)
mate2 <- read.table(paste(path2data,"/s_2_sequence.mate2.avgQualPerBP.txt",sep=""),header=TRUE)
     
# mate 1

jpeg("s_2_1_sequence.avgQualPerBP.jpeg")

plot(mate1[,1],mate1[,2],xlab="Base pair",ylab="Average quality",main="Mate 1",pch=20,ylim=c(-5,40))
arrows(mate1[,1],mate1[,2]-mate1[,3],mate1[,1],mate1[,2]+mate1[,3],length=0.05,angle=90,code=3,lty=1,lwd=0.3)
abline(h=20,lty=3)

dev.off()

# mate 2

jpeg("s_2_2_sequence.avgQualPerBP.jpeg")

plot(mate2[,1],mate2[,2],xlab="Base pair",ylab="",main="Mate 2",pch=20,ylim=c(-5,40))
arrows(mate2[,1],mate2[,2]-mate2[,3],mate2[,1],mate2[,2]+mate2[,3],length=0.05,angle=90,code=3,lty=1,lwd=0.3)
abline(h=20,lty=3)

dev.off()
