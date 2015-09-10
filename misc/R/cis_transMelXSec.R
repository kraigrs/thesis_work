#######################################################################
#
# Name: Binom-FET-PoO_070510-try1.R
#
# Author: Written by Joe Coolon, Kraig Stevenson and Joel McManus 
#
# Date: 10/19/09 (modified 070510, and 7/06/10, 07/07/2010)
#       08/10/2010: modified to only take 1 cross direction
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

data <- read.delim(commandArgs()[3],header=TRUE,sep="\t");
#data <- read.delim("/Users/kraigrs/Wittkopp/mel_mel_data/zhr_z30.counts.txt",header=TRUE,sep="\t");

# Run Binomial Exact Tests on the hybrids using the Pearson-Clopper Method 
# to determine cis-regulatory changes. Change output file names to whatever you want. 

# initiate variables

pvalsHyb <- NULL; pvalsHyb_Par <- NULL;
pvalsPar <- NULL; 

tempH <- NULL; tempHw <- NULL; tempH_P <- NULL;
tempP <- NULL; tempPw <- NULL;

for (i in 1:nrow(data)) 
{
	Hyb_Bi <- binom.test(round(data$hybS1[[i]]), round(data$hybTot[[i]]), p = 0.5, alternative = c("t"), conf.level = 0.95);
	Par_Bi <- binom.test(round(data$parS1[[i]]), round(data$parTot[[i]]), p = 0.5, alternative = c("t"), conf.level = 0.95);
	
	# collect p-values from binomial tests
	pvalsHyb <- rbind(pvalsHyb,Hyb_Bi$p.value);
	pvalsPar <- rbind(pvalsPar,Par_Bi$p.value);
	
	tempH <- rbind(tempH,c(Hyb_Bi$estimate,Hyb_Bi$conf.int[1],Hyb_Bi$conf.int[2],Hyb_Bi$p.value));
	tempP <- rbind(tempP,c(Par_Bi$estimate,Par_Bi$conf.int[1],Par_Bi$conf.int[2],Par_Bi$p.value));
	
	# collect Wilson confidence intervals for binomial tests (to be used later for plotting seq. vs pyro)
	Hyb_Bi2 <- binconf(round(data$hybS1[[i]]), round(data$hybTot[[i]]), alpha = 0.05, method = c("wilson"));
	Par_Bi2 <- binconf(round(data$parS1[[i]]), round(data$parTot[[i]]), alpha = 0.05, method = c("wilson"));
	
	tempHw <- rbind(tempHw,c(Hyb_Bi2[2],Hyb_Bi2[3]));
	tempPw <- rbind(tempPw,c(Par_Bi2[2],Par_Bi2[3]));
	
	# Fisher's exact test from 2x2 tables of counts
	hyb_par.FET <- fisher.test(matrix(c(round(data$hybS1[[i]]),round(data$hybS2[[i]]),round(data$parS1[[i]]),round(data$parS2[[i]])),nr=2));

	# collect p-values from FETs
	pvalsHyb_Par <- rbind(pvalsHyb_Par,hyb_par.FET$p.value);
	
	tempH_P <- rbind(tempH_P,c(hyb_par.FET$estimate,hyb_par.FET$conf.int[1],hyb_par.FET$conf.int[2],hyb_par.FET$p.value));
}

pvalsHyb.adj <- p.adjust(pvalsHyb,method="fdr");
pvalsPar.adj <- p.adjust(pvalsPar,method="fdr");
	
pvalsHyb_Par.adj <- p.adjust(pvalsHyb_Par,method="fdr");
               
data.output <- cbind(data[,1:4],tempH,pvalsHyb.adj,tempHw,data[,5:7],tempP,pvalsPar.adj,tempPw,tempH_P,pvalsHyb_Par.adj);

colnames(data.output) <- c("gene","H.s1","H.s2","H.tot","H.binomEst","H.LB","H.UB","H.binomP","H.binomQ","H.LBw","H.UBw","P.s1","P.s2","P.tot","P.binomEst","P.LB","P.UB","P.binomP","P.binomQ","P.LBw","P.UBw","H_P.fetEst","H_P.LB","H_P.UB","H_P.P","H_P.Q");
	
write.table(data.output,file=commandArgs()[4],sep="\t",quote=FALSE,row.names=FALSE);
