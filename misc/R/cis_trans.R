#######################################################################
#
# Name: Binom-FET-PoO_070510-try1.R
#
# Author: Written by Joe Coolon, Kraig Stevenson and Joel McManus 
#
# Date: 10/19/09 (modified 070510, and 7/06/10, 07/07/2010)
#
# Purpose: This code will be used to first perform Binomial Exact 
# Tests using the Wilson method.  Using this test we will test
# for differential expression between parent species as well
# as test for differences in allele-specific expression in 
# hybrids.  Additionally this script will perform Fisher's Exact
# Tests comparing parental and hybrid expression. Fianlly, this script 
# will perform Fisher's Eact tests to determoine parent of origin effects.
#
# Details: 
#
# Input: Count data from perl script XXXX.pl that has counted 
# of times a particular allele (snp) was identified in sequencing 
# output.
#
#
# Output: File containing p-values, confidence interval and estimate.
#
#
# Usage: Will be used following scripts that process Illumina 
# sequencing data into allele-specific counts.
#
#######################################################################

# load packages to be used.

library(binom)
library(Hmisc)

# load data files to be used, change below according to the names of 
# files intended to be used. For parent of origin tests need two hyb files.

#data <- read.delim("/Users/wittkopp-lab/Desktop/Kraig/Wittkopp/mockdata/data.txt");
data <- read.delim(commandArgs()[3],header=TRUE,sep="\t");

# Run Binomial Exact Tests on the hybrids using the Pearson-Clopper Method 
# to determine cis-regulatory changes. Change output file names to whatever you want. 
# Tests for both directions of cross are included (hyb1 and hyb2).

# initiate variables

pvalsHyb1 <- NULL; pvalsHyb1_Par <- NULL;
pvalsHyb2 <- NULL; pvalsHyb2_Par <- NULL;
pvalsPar <- NULL; pvalsHyb1_Hyb2 <- NULL;

tempH1 <- NULL; tempH1w <- NULL; tempH1_P <- NULL;
tempH2 <- NULL; tempH2w <- NULL; tempH2_P <- NULL;
tempP <- NULL; tempPw <- NULL; tempH1_H2 <- NULL;

for (i in 1:nrow(data)) 
{
	Hyb1_Bi <- binom.test(round(data$hyb1s1[[i]]), round(data$hyb1tot[[i]]), p = 0.5, alternative = c("t"), conf.level = 0.95);
	Hyb2_Bi <- binom.test(round(data$hyb2s1[[i]]), round(data$hyb2tot[[i]]), p = 0.5, alternative = c("t"), conf.level = 0.95);
	Par_Bi <- binom.test(round(data$parS1[[i]]), round(data$parTot[[i]]), p = 0.5, alternative = c("t"), conf.level = 0.95);
	
	# collect p-values from binomial tests
	pvalsHyb1 <- rbind(pvalsHyb1,Hyb1_Bi$p.value);
	pvalsHyb2 <- rbind(pvalsHyb2,Hyb2_Bi$p.value); 
	pvalsPar <- rbind(pvalsPar,Par_Bi$p.value);
	
	tempH1 <- rbind(tempH1,c(Hyb1_Bi$estimate,Hyb1_Bi$conf.int[1],Hyb1_Bi$conf.int[2],Hyb1_Bi$p.value));
	tempH2 <- rbind(tempH2,c(Hyb2_Bi$estimate,Hyb2_Bi$conf.int[1],Hyb2_Bi$conf.int[2],Hyb2_Bi$p.value));
	tempP <- rbind(tempP,c(Par_Bi$estimate,Par_Bi$conf.int[1],Par_Bi$conf.int[2],Par_Bi$p.value));
	
	# collect Wilson confidence intervals for binomial tests (to be used later for plotting seq. vs pyro)
	Hyb1_Bi2 <- binconf(round(data$hyb1s1[[i]]), round(data$hyb1tot[[i]]), alpha = 0.05, method = c("wilson"));
	Hyb2_Bi2 <- binconf(round(data$hyb2s1[[i]]), round(data$hyb2tot[[i]]), alpha = 0.05, method = c("wilson"));
	Par_Bi2 <- binconf(round(data$parS1[[i]]), round(data$parTot[[i]]), alpha =0.05, method = c("wilson"));
	
	tempH1w <- rbind(tempH1w,c(Hyb1_Bi2[2],Hyb1_Bi2[3]));
	tempH2w <- rbind(tempH2w,c(Hyb2_Bi2[2],Hyb2_Bi2[3]));
	tempPw <- rbind(tempPw,c(Par_Bi2[2],Par_Bi2[3]));
	
	# Fisher's exact test from 2x2 tables of counts
	hyb1_par.FET <- fisher.test(matrix(c(round(data$hyb1s1[[i]]),round(data$hyb1s2[[i]]),round(data$parS1[[i]]),round(data$parS2[[i]])),nr=2));
	hyb2_par.FET <- fisher.test(matrix(c(round(data$hyb2s1[[i]]),round(data$hyb2s2[[i]]),round(data$parS1[[i]]),round(data$parS2[[i]])),nr=2));
	hyb1_hyb2.FET <- fisher.test(matrix(c(round(data$hyb1s1[[i]]),round(data$hyb1s2[[i]]),round(data$hyb2s1[[i]]),round(data$hyb2s2[[i]])),nr=2));
	
	# collect p-values from FETs
	pvalsHyb1_Par <- rbind(pvalsHyb1_Par,hyb1_par.FET$p.value);
	pvalsHyb2_Par <- rbind(pvalsHyb2_Par,hyb2_par.FET$p.value); 
	pvalsHyb1_Hyb2 <- rbind(pvalsHyb1_Hyb2,hyb1_hyb2.FET$p.value);
	
	tempH1_P <- rbind(tempH1_P,c(hyb1_par.FET$estimate,hyb1_par.FET$conf.int[1],hyb1_par.FET$conf.int[2],hyb1_par.FET$p.value));
	tempH2_P <- rbind(tempH2_P,c(hyb2_par.FET$estimate,hyb2_par.FET$conf.int[1],hyb2_par.FET$conf.int[2],hyb2_par.FET$p.value));
	tempH1_H2 <- rbind(tempH1_H2,c(hyb1_hyb2.FET$estimate,hyb1_hyb2.FET$conf.int[1],hyb1_hyb2.FET$conf.int[2],hyb1_hyb2.FET$p.value));
}

pvalsHyb1.adj <- p.adjust(pvalsHyb1,method="fdr");
pvalsHyb2.adj <- p.adjust(pvalsHyb2,method="fdr"); 
pvalsPar.adj <- p.adjust(pvalsPar,method="fdr");
	
pvalsHyb1_Par.adj <- p.adjust(pvalsHyb1_Par,method="fdr");
pvalsHyb2_Par.adj <- p.adjust(pvalsHyb2_Par,method="fdr"); 
pvalsHyb1_Hyb2.adj <- p.adjust(pvalsHyb1_Hyb2,method="fdr");

data.output <- cbind(data[,1:4],tempH1,pvalsHyb1.adj,tempH1w,data[,5:7],tempH2, pvalsHyb2.adj,tempH2w,data[,8:10],tempP,pvalsPar.adj,tempPw,tempH1_P,pvalsHyb1_Par.adj,tempH2_P,pvalsHyb2_Par.adj,tempH1_H2,pvalsHyb1_Hyb2.adj);

colnames(data.output) <- c("gene","H1.s1","H1.s2","H1.tot","H1.binomEst","H1.LB","H1.UB","H1.binomP","H1.binomQ","H1.LBw","H1.UBw","H2.s1","H2.s2","H2.tot","H2.binomEst","H2.LB","H2.UB","H2.binomP","H2.binomQ","H2.LBw","H2.UBw","P.s1","P.s2","P.tot","P.binomEst","P.LB","P.UB","P.binomP","P.binomQ","P.LBw","P.UBw","H1_P.fetEst","H1_P.LB","H1_P.UB","H1_P.P","H1_P.Q","H2_P.fetEst","H2_P.LB","H2_P.UB","H2_P.P","H2_P.Q","H1_H2.fetEst","H1_H2.LB","H1_H2.UB","H1_H2.P","H1_H2.Q");
	
write.table(data.output,file=commandArgs()[4],sep="\t",quote=FALSE,row.names=FALSE);