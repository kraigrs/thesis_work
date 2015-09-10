######################################################################
#
# Name: Dominance.R
#
# Author: Wirtten by Trisha Wittkopp, Joel McManus and modified 
# by Joe Coolon. Based on analysis by Greg Gibson (Gibson et al. 2004)
#
# Date: 07/05/10 (modified, 10/08/09 originally written)
#
# Purpose: This code will be used to evaluate the mode of 
# inheretance of the total experssion levels (counts) from
# kraigware X.X output.  The Dominance classifications follow 
# those from Gibson et al. 2004 
#
# Details: "Total expression is normalized by dividing the number 
# of mapped reads at each gene by the total number of mapped reads 
# (needs to be altered to apply to the particular sequencing results 
# being analyzed) for the entire genome, and multiplying by 100 to 
# give a percent expression value. Log-transformed percent expression 
# values of each parent are then subtracted from those of the hybrid 
# to examine changes in expression. Genes whose total expression in 
# hybrids deviates more than 1.25-fold from that of either parent are 
# considered to have non-conserved inheritance, and are classified  
# as having additive, dominant, underdominant, or overdominant  
# inheritance based on the magnitude of the difference between  
# total expression in the hybrid and in each parental species.
# Briefly, genes for which expression in the hybrid was less 
# than Parent species 1 and greater than Parent species 2 (or vice 
# versa) were classified as additive; genes for which expression 
# in the hybrid was similar to one of the parents were classified 
# as dominant; and genes for which expression in the hybrid was 
# either greater than or less than both Parent species 1 and 
# Parent species 2 were considered misexpressed and classified 
# as over-dominant and under-dominant, respectively. Regardless 
# of statistical significance, any pair of genotypes with 
# expression that differed less than 1.25-fold was considered 
# to have similar expression for this analysis (Gibson et al. 2004)."
# (Mcmanus et al 2010).
#
# Input: Count data from perl script XXXX.pl that has counted 
# of times a particular allele (snp) was identified in sequencing 
# output from both mixed parentals and hybrids.
#
#
# Output: Plot of mode of inheretance that is color-coded. 
#
#
# Usage: Will be used following scripts that process Illumina 
# sequencing data into allele-specific counts.
#
########################################################################

input <- commandArgs()[3];

# Read in the data (change below to reflect the datafile you want to use
data <- read.delim(input,header=TRUE);

# Restrict to genes that are significantly differentiallly expressed 
# between the two parental species.

div <- subset(data, "variable of mixed q values" < 0.05);
attach(div);

# Convert to percentage of total
Hyb_pct <- Hyb_Mel.Sec/(33358462);

# Make plot of parental vs hybrid
plot(log2(Hyb_pct), log2(Mix_M_pct));
points(log2(Hyb_pct),log2(Mix_S_pct), col="red3");

# Calculate values for axis of dominance plot
mdiff <- log2(Mix_M_pct) - log2(Hyb_pct*100);
sdiff <- log2(Mix_S_pct) - log2(Hyb_pct*100);

# Make Dominance plot
plot(mdiff,sdiff);

# Perform tests for determining classifications for each gene
test <- cbind(mdiff,sdiff);
a <- subset(test,test[,1] < -0.097 & test[,2] < -0.097);
b <- subset(test,test[,1] < -0.097 & test[,2] > -0.097 & test[,2] <  0.097);
c <- subset(test,test[,1] < -0.097 & test[,2] >  0.097);
d <- subset(test,test[,1] > -0.097 & test[,1] <  0.097 & test[,2] < -0.097);
e <- subset(test,test[,1] > -0.097 & test[,1] <  0.097 & test[,2] > -0.097 & test[,2] <  0.097);
f <- subset(test,test[,1] > -0.097 & test[,1] <  0.097 & test[,2] >  0.097);
g <- subset(test,test[,1] >  0.097 & test[,2] < -0.097);
h <- subset(test,test[,1] >  0.097 & test[,2] > -0.097 & test[,2] <  0.097);
i <- subset(test,test[,1] >  0.097 & test[,2] >  0.097);
temp <- c(nrow(a),nrow(b),nrow(c),nrow(d),nrow(e),nrow(f),nrow(g),nrow(h),nrow(i));
cnts <- matrix(data=temp,nrow=3, ncol=3,byrow=FALSE, dimnames=list(c("<", "=(<1.25 fold)", ">"), c("<", "=(<1.25 fold)", ">")));
#cnts
plot(test, type="n", xlab="log2(Mel exp) - log2(Hyb exp)", ylab="log2(Sec exp) - log2(Hyb exp)");
#overdominant
points(a, col="4");
#additive
points(c, col="5");
points(g, col="5");
#underdominant
points(i, col="2");
#mel dominant
points(d, col="3");
points(f, col="3");
#sim dominant
points(b ,col="6");
points(h, col="6");
#conserved
points(e, col="orange");
#nrow(data)-nrow(div);