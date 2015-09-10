################################################################################################################
# 
# 11/16/2009
# 
# histReads.R
#
# Purpose: input a .txt file from readInfo.pl containing the best gene hit, exon number, and counts for
#          melanogaster and simulans and contruct histograms
#
# Output: a histogram showing the log2(M/S) ration of mel to sim reads
#
# Example: R --slave < histReads.R infile outfile_suffix > ./trash
#
################################################################################################################

data <- read.delim(commandArgs()[3],header=TRUE,sep="\t")
#data <- read.delim("filename",header=TRUE,sep="\t")

file_suffix <- commandArgs()[4]

mel <- data[,5]
sim <- data[,6]

log2ratio <- log2(mel/sim)
log2ratio <- log2ratio[!(is.na(log2ratio)|is.infinite(log2ratio))]

jpeg(paste(file_suffix,".melSimRatio.jpeg",sep=""))
hist(log2ratio,main="Distribution of log2(M/S)",xlab="log2(M/S)",right=FALSE)
legend("topright",legend=c(paste("Mean: ",round(mean(log2ratio),digits=2),sep=""),
                           paste("95% CI: [",round(quantile(log2ratio,.025),digits=2),",",
                           round(quantile(log2ratio,.975),digits=2),"]",sep="")))
dev.off()
