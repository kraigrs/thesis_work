################################################################################################################
# 
# 10/15/2009
# 
# histFASTQ.R
#
# Purpose: input the output from summstatFASTQ in order to make some nice
#          histograms, show correlation, and linear regression
#
# Output: 2 histograms showing the distribution of the maximum
#         cluster enrichment scores and the averages of the cluster maxima
#
# Example: R --slave < histFASTQ.R s_2_1_sequence.summary.txt s_2_2_sequence.summary.txt s_2_sequence > ./trash
#
################################################################################################################

mate1 <- read.delim(commandArgs()[3],header=TRUE,sep="\t")
#mate1 <- read.delim("s_2_1_sequence.summary.txt",header=TRUE,sep="\t")
mate2 <- read.delim(commandArgs()[4],header=TRUE,sep="\t")
#mate2 <- read.delim("s_2_2_sequence.summary.txt",header=TRUE,sep="\t")

lane <- commandArgs()[5]
#lane <- "s_2_sequence"

n1 <- mate1[,2]
q1 <- mate1[,3]

n2 <- mate2[,2]
q2 <- mate2[,3]

jpeg(paste(lane,"Nmate1.jpeg",sep=""))
hist(n1,main="Distribution of ambiguous bases per read for mate1",xlab="N count",right=FALSE)
legend("topright",legend=c(paste("Mean: ",round(mean(n1),digits=2),sep=""),
                           paste("95% CI: [",round(quantile(n1,.025),digits=2),",",
                           round(quantile(n1,.975),digits=2),"]",sep="")))
dev.off()

jpeg(paste(lane,"Nmate2.jpeg",sep=""))
hist(n2,main="Distribution of ambiguous bases per read for mate2",xlab="N count",right=FALSE)
legend("topright",legend=c(paste("Mean: ",round(mean(n2),digits=2),sep=""),
                           paste("95% CI: [",round(quantile(n2,.025),digits=2),",",
                           round(quantile(n2,.975),digits=2),"]",sep="")))
dev.off()

jpeg(paste(lane,"Qmate1.jpeg",sep=""))
hist(q1,main="Distribution of average read quality mate1",xlab="Avg. quality",right=FALSE)
legend("topright",legend=c(paste("Mean: ",round(mean(q1),digits=2),sep=""),
                           paste("95% CI: [",round(quantile(q1,.025),digits=2),",",
                           round(quantile(q1,.975),digits=2),"]",sep="")))
dev.off()

jpeg(paste(lane,"Qmate2.jpeg",sep=""))
hist(q2,main="Distribution of average read quality mate2",xlab="Avg. quality",right=FALSE)
legend("topright",legend=c(paste("Mean: ",round(mean(q2),digits=2),sep=""),
                           paste("95% CI: [",round(quantile(q2,.025),digits=2),",",
                           round(quantile(q2,.975),digits=2),"]",sep="")))
dev.off()

jpeg(paste(lane,"QvsNmate1.jpeg",sep=""))
plot(n1,q1,main="Average quality score vs. ambiguous base count for mate 1")
reg <- lm(q1~n1)
abline(a=reg$coefficients[1],b=reg$coefficients[2])
dev.off()

jpeg(paste(lane,"QvsNmate2.jpeg",sep=""))
plot(n2,q2,main="Average quality score vs. ambiguous base count for mate 2")
reg <- lm(q2~n2)
abline(a=reg$coefficients[1],b=reg$coefficients[2])
dev.off()

jpeg(paste(lane,"N2vsN1.jpeg",sep=""))
plot(n1,n2,main="Correlation of ambiguous bases for each mate")
reg <- lm(n2~n1)
abline(a=0,b=1)
dev.off()

jpeg(paste(lane,"Q2vsQ1.jpeg",sep=""))
plot(q1,q2,main="Correlation of average quality scores for each mate")
reg <- lm(q2~q1)
abline(a=0,b=1)
dev.off()











