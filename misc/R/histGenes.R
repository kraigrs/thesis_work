################################################################################################################
# 
# 11/17/2009
# 
# histGenes.R
#
# Purpose: input a .txt file from geneSummary.pl containing the number of SNPs, reads, and M vs S counts on each gene
#
# Output: histograms showing the distributions of those things
#
# Example: R --slave < histGenes.R infile outfile_suffix > ./trash
#
################################################################################################################

data <- read.delim(commandArgs()[3],header=TRUE,sep="\t")
#data <- read.delim("filename",header=TRUE,sep="\t")

file_suffix <- commandArgs()[4]

data <- data[!(data[,2]==0),]

SNPs <- data[,2]
reads <- data[,3]

jpeg(paste(file_suffix,".SNPsHisto.jpeg",sep=""))
hist(SNPs,main="Distribution of SNPs per gene",xlab="Number of SNPs per gene",right=FALSE)
legend("topright",legend=c(paste("Mean: ",round(mean(SNPs),digits=2),sep=""),
                           paste("95% CI: [",round(quantile(SNPs,.025),digits=2),",",
                           round(quantile(SNPs,.975),digits=2),"]",sep="")))
dev.off()

jpeg(paste(file_suffix,".readsHisto.jpeg",sep=""))
hist(reads,main="Distribution of reads per gene",xlab="Number of reads per gene",right=FALSE)
legend("topright",legend=c(paste("Mean: ",round(mean(reads),digits=2),sep=""),
                           paste("95% CI: [",round(quantile(reads,.025),digits=2),",",
                           round(quantile(reads,.975),digits=2),"]",sep="")))
dev.off()

