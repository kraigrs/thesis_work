######################################################################
#
# Name: cis-trans_plot.R
#
# Author: Wirtten by Joe Coolon. 
#
# Date: 08/03/10
#
# Purpose: 
#
# Details: 
#
# Input: Count data from perl script classify.pl that has counted 
# of times a particular allele (snp) was identified in sequencing 
# output from both mixed parentals and hybrids.
#
#
# Output: Plot of hybrid vs parental that is color-coded. 
#
#
# Usage: Will be used following scripts that process Illumina 
# sequencing data into allele-specific counts.
#
########################################################################


# Read in the data (change below to reflect the datafile you want to use
#data <- read.delim(commandArgs()[3],header=TRUE,sep="\t");
data <- read.delim("/Users/kraigrs/Wittkopp/mel_sec_data/mel_sec.divergence.txt",header=TRUE,sep="\t");

# Make zhr/z30 ratios fopr both parental and hybrid
hyb_ratio <- (data$hybS1/data$hybS2);
Par_ratio <- (data$parS1/data$parS2);

#Log2 transfor the ratios
log_hyb_ratio <- log2(hyb_ratio);
log_Par_ratio <- log2(Par_ratio);

ratios <- NULL;

for(i in 1:nrow(data))
{
	temp <- cbind(log_Par_ratio[i], log_hyb_ratio[i], toString(data$class[i]));
	ratios <- rbind(ratios, temp);
}

# Perform tests for determining classifications for each gene

cis <- subset(ratios,ratios[,3] == "cis");
trans <- subset(ratios,ratios[,3] == "trans");
cisPtrans <- subset(ratios,ratios[,3] == "cis+trans");
cisXtrans <- subset(ratios,ratios[,3] == "cisXtrans");
comp <- subset(ratios,ratios[,3] == "compensatory");
cons <- subset(ratios,ratios[,3] == "conserved");
amb <- subset(ratios,ratios[,3] == "ambiguous");

plot(log_Par_ratio, log_hyb_ratio, type="n", xlab="log2(parental(mel/sec))", ylab="log2(hybrid(mel/sec))",main="melXsec",xlim=c(-10,10),ylim=c(-10,10))
#plot(log_Par_ratio, log_hyb_ratio, type="n", xlab="log2(parental(mel/sec))", ylab="log2(hybrid(mel/sec))",xlim=c(-2,2),ylim=c(-2,2))
#plot(log_Par_ratio, log_hyb_ratio, type="p", xlab="log2(parental(mel/sec))", ylab="log2(hybrid(mel/sec))")

#Ambiguous
points(amb, col=rgb(t(col2rgb("gray")),alpha=50,maxColorValue=255),pch=19,cex=0.25)

#cis + trans
points(cisPtrans, col=rgb(t(col2rgb("blue")),alpha=50,maxColorValue=255),pch=19,cex=0.25)

#cis X trans
points(cisXtrans, col=rgb(t(col2rgb("green")),alpha=50,maxColorValue=255),pch=19,cex=0.25)

#All cis
points(cis, col=rgb(t(col2rgb("black")),alpha=50,maxColorValue=255),pch=19,cex=0.25)

#All trans
points(trans, col=rgb(t(col2rgb("red")),alpha=50,maxColorValue=255),pch=19,cex=0.25)

#Compensatory
points(comp, col=rgb(t(col2rgb("orange")),alpha=50,maxColorValue=255),pch=19,cex=0.25)

#Conserved
points(cons, col=rgb(t(col2rgb("yellow")),alpha=50,maxColorValue=255),pch=19,cex=0.25)

abline(h=0,lty=2)
abline(v=0,lty=2)
#abline(a=0,b=1,lty=2)
