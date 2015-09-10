############################################################################################
# script to characterize regulatory divergence
############################################################################################

#Usage: R --vanilla --args parent hybrid dir out < regulatory_divergence.R > regulatory_divergence.out

############################################################################################
# functions

# test for parental and hybrid differences in allele1 and allele2
FET_P <- function(locus)
{
	table <- matrix(c(locus[1],locus[3]-locus[1],locus[2],locus[4]-locus[2]),
		            nrow=2,
		            dimnames=list(c("Gene","Not gene"),c("Par1","Par2")));
	test <- fisher.test(table,or=1,alternative="two.sided",conf.level=0.95);
	result <- c(test$estimate,test$conf.int[1],test$conf.int[2],test$p.value);
	return(result);
}

# test for hybrid differences in allele1 and allele2
BET <- function(locus)
{
	test <- binom.test(locus[1],locus[1]+locus[2],0.5,alternative="two.sided",conf.level=0.95);
	result <- c(test$estimate,test$conf.int[1],test$conf.int[2],test$p.value);
	return(result);
}

# test for parental and hybrid differences in allele1 and allele2
FET_PH <- function(locus)
{
	table <- matrix(c(locus[1],locus[3],locus[2],locus[4]),
					nrow=2,
					dimnames=list(c("Parents","Hybrid"),c("Allele1","Allele2")));
	test <- fisher.test(table,or=1,alternative="two.sided",conf.level=0.95);
	result <- c(test$estimate,test$conf.int[1],test$conf.int[2],test$p.value);
	return(result);
}

# regulatory divergence classification
regulatory_divergence <- function(locus)
{
	p1 <- locus[1];
	p2 <- locus[2];
	h1 <- locus[3];
	h2 <- locus[4];

	pdiff <- locus[5];
	hdiff <- locus[6];
	phdiff <- locus[7];
	
	if(pdiff >= 0.05 & hdiff >= 0.05 & phdiff >= 0.05){return("conserved");}
	else if(pdiff < 0.05 & hdiff < 0.05 & phdiff >= 0.05){return("all_cis");}
	else if(pdiff < 0.05 & hdiff >= 0.05 & phdiff < 0.05){return("all_trans");}
	else if(pdiff < 0.05 & hdiff < 0.05 & phdiff < 0.05)
	{
		if(log2(p1/p2)/log2(h1/h2) > 1){return("cisPtrans");}
		else if(log2(p1/p2)/log2(h1/h2) < 1){return("cisXtrans");}
	}
	else if(pdiff >= 0.05 & hdiff < 0.05 & phdiff < 0.05){return("compensatory");}
	else if(pdiff < 0.05 & hdiff >= 0.05 & phdiff >= 0.05){return("ambiguous");}
	else if(pdiff >= 0.05 & hdiff < 0.05 & phdiff >= 0.05){return("ambiguous");}
	else if(pdiff >= 0.05 & hdiff >= 0.05 & phdiff < 0.05){return("ambiguous");}
}

############################################################################################
# read in the data, apply cutoff(s) and filters, perform functions

args <- commandArgs(trailingOnly = TRUE);
#print(args);
data <- read.table(args[1],header=TRUE,sep="\t");
dir <- args[2];
out <- args[3];

#data <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_F_carcass.genes_down_ASE.txt",header=TRUE,sep="\t");
#dir <- "/Users/kraigrs/Wittkopp/Machado_Dpse/plots";
#out <- "TL_Toro1_H6_F_carcass";

cutoff <- 20;

reliable <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/reliable_genes.txt",header=TRUE,sep="\t");
temp <- merge(data,reliable,by.x="gene",by.y="gene");

data_cutoff <- temp[which(temp[2]+temp[3] >= cutoff & temp[4]+temp[5] >= cutoff),];


# parental counts for species1 and species2
#FET_P_results <- t(apply(as.matrix(cbind(data_cutoff[,c(2,3)],
#						   rep(sum(data_cutoff[,2]),nrow(data_cutoff)),
#						   rep(sum(data_cutoff[,3]),nrow(data_cutoff)))),1,FET_P));
#colnames(FET_P_results) <- c("FET_P_estimate","FET_P_lower","FET_P_upper","FET_P_pval");

BET_P_results <- t(apply(as.matrix(data_cutoff[,c(2,3)]),1,BET));				
colnames(BET_P_results) <- c("BET_P_estimate","BET_P_lower","BET_P_upper","BET_P_pval");

# hybrid counts for species1 and species2
BET_H_results <- t(apply(as.matrix(data_cutoff[,c(4,5)]),1,BET));
colnames(BET_H_results) <- c("BET_H_estimate","BET_H_lower","BET_H_upper","BET_H_pval");

# parental counts for species1 and species2, then hybrid counts for species1 and species2
FET_PH_results <- t(apply(as.matrix(data_cutoff[,c(2,3,4,5)]),1,FET_PH));
colnames(FET_PH_results) <- c("FET_PH_estimate","FET_PH_lower","FET_PH_upper","FET_PH_pval");

data_results <- cbind(data_cutoff,BET_P_results,BET_H_results,FET_PH_results);

BET_P_qval <- p.adjust(data_results$BET_P_pval,method="fdr");
BET_H_qval <- p.adjust(data_results$BET_H_pval,method="fdr");
FET_PH_qval <- p.adjust(data_results$FET_PH_pval,method="fdr");

data_results <- cbind(data_results,BET_P_qval,BET_H_qval,FET_PH_qval);

temp <- data_results[!is.na(log2(data_results[,2]/data_results[,3])/log2(data_results[,4]/data_results[,5])),];
data_results <- temp;

percent_cis <- abs(log2(data_results[,4]/data_results[,5]))/(abs(log2(data_results[,4]/data_results[,5]))+abs(log2(data_results[,2]/data_results[,3])-log2(data_results[,4]/data_results[,5])));

reg_div <- apply(as.matrix(data_results[,c(2,3,4,5,18,19,20)]),1,regulatory_divergence);

data_final <- cbind(data_results,percent_cis,reg_div);

write.table(data_final,file=paste(dir,"/",out,"_reg_div.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE);

# violin plot of percent cis by regulatory divergence category
library(vioplot);

cis <- data_final$percent_cis[data_final$reg_div=="all_cis"];
trans <- data_final$percent_cis[data_final$reg_div=="all_trans"];
cisPtrans <- data_final$percent_cis[data_final$reg_div=="cisPtrans"];
cisXtrans <- data_final$percent_cis[data_final$reg_div=="cisXtrans"];
compensatory <- data_final$percent_cis[data_final$reg_div=="compensatory"];
conserved <- data_final$percent_cis[data_final$reg_div=="conserved"];
ambiguous <- data_final$percent_cis[data_final$reg_div=="ambiguous"];

data_final$reg_div <- factor(data_final$reg_div,levels(data_final$reg_div)[c(1,2,4,5,6,7,3)]);

avg <- mean(data_final$percent_cis[which(is.finite(data_final$percent_cis))]);

temp <- table(data_final$reg_div);
parms <- scale(sqrt(temp),center=FALSE);
cats <- names(temp);

pdf(file=paste(dir,"/",out,"_percent_cis_reg_div.pdf",sep=""));

plot(0:12,rep(avg,13),type='l',xlim=c(0,12),ylim=c(0,1),ylab="",xlab="",xaxt="n",lty=3,main=paste(out," %cis",sep=""));

if(length(cis) > 0){vioplot(cis[is.finite(cis)],at=1,wex=parms[which(cats=="all_cis")],add=TRUE,names="",col="gray20",pchMed=19,colMed="gray40");}
if(length(trans) > 0){vioplot(trans[is.finite(trans)],at=2.5,wex=parms[which(cats=="all_trans")],add=TRUE,names="",col="red",pchMed=19,colMed="gray40");}
if(length(cisPtrans) > 0){vioplot(cisPtrans[is.finite(cisPtrans)],at=4,wex=parms[which(cats=="cisPtrans")],add=TRUE,names="",col="purple",pchMed=19,colMed="gray40");}
if(length(cisXtrans) > 0){vioplot(cisXtrans[is.finite(cisXtrans)],at=5.5,wex=parms[which(cats=="cisXtrans")],add=TRUE,names="",col="green",pchMed=19,colMed="gray40");}
if(length(compensatory) > 0){vioplot(compensatory[is.finite(compensatory)],at=7,wex=parms[which(cats=="compensatory")],add=TRUE,names="",col="orange",pchMed=19,colMed="gray40");}
if(length(conserved) > 0){vioplot(conserved[is.finite(conserved)],at=8.5,wex=parms[which(cats=="conserved")],add=TRUE,names="",col="yellow",pchMed=19,colMed="gray40");}
if(length(ambiguous) > 0){vioplot(ambiguous[is.finite(ambiguous)],at=10,wex=parms[which(cats=="ambiguous")],add=TRUE,names="",col="gray",pchMed=19,colMed="gray40");}
vioplot(data_final$percent_cis[is.finite(data_final$percent_cis)],at=11.5,wex=parms[which(cats=="ambiguous")],add=TRUE,names="",col="white",pchMed=19,colMed="gray40");

dev.off();

data_rm <- subset(data_final,data_final[,2]>0&data_final[,3]>0&data_final[,4]>0&data_final[,5]>0);

data_rm$reg_div <- factor(data_rm$reg_div,levels=c("all_cis","all_trans","cisPtrans","cisXtrans","compensatory","conserved","ambiguous"));

cats <- table(data_rm$reg_div);

pdf(file=paste(dir,"/",out,"_reg_div_pie.pdf",sep=""));
pie(c(cats[1],cats[3]+cats[4]+cats[5],cats[2],cats[6],cats[7]),labels="",col=c("blue","purple","red","black","gray"),
    clockwise=TRUE,init.angle=180);
dev.off();

#axis <- ceiling(max(c(max(abs(log2(data_rm[,2]/data_rm[,3]))),max(abs(log2(data_rm[,4]/data_rm[,5]))))));
axis <- 5;

cis <- data_rm[data_rm$reg_div=="all_cis",];
trans <- data_rm[data_rm$reg_div=="all_trans",];
cisPtrans <- data_rm[data_rm$reg_div=="cisPtrans",];
cisXtrans <- data_rm[data_rm$reg_div=="cisXtrans",];
compensatory <- data_rm[data_rm$reg_div=="compensatory",];
conserved <- data_rm[data_rm$reg_div=="conserved",];
ambiguous <- data_rm[data_rm$reg_div=="ambiguous",];

pdf(file=paste(dir,"/",out,"_reg_div.pdf",sep=""));

plot(log2(conserved[,2]/conserved[,3]),log2(conserved[,4]/conserved[,5]),xlim=c(-axis,axis),ylim=c(-axis,axis),col=rgb(t(col2rgb("black")),alpha=100,maxColorValue=255),pch=19,cex=0.25,xlab="log2(A1/A2) parents",ylab="log2(A1/A2) hybrid",main=paste(out,sep=""));

points(log2(ambiguous[,2]/ambiguous[,3]),log2(ambiguous[,4]/ambiguous[,5]),xlim=c(-axis,axis),ylim=c(-axis,axis),col=rgb(t(col2rgb("gray")),alpha=100,maxColorValue=255),pch=19,cex=0.25);

points(log2(compensatory[,2]/compensatory[,3]),log2(compensatory[,4]/compensatory[,5]),xlim=c(-axis,axis),ylim=c(-axis,axis),col=rgb(t(col2rgb("purple")),alpha=100,maxColorValue=255),pch=19,cex=0.25);

points(log2(cisXtrans[,2]/cisXtrans[,3]),log2(cisXtrans[,4]/cisXtrans[,5]),xlim=c(-axis,axis),ylim=c(-axis,axis),col=rgb(t(col2rgb("purple")),alpha=100,maxColorValue=255),pch=19,cex=0.25);

points(log2(cisPtrans[,2]/cisPtrans[,3]),log2(cisPtrans[,4]/cisPtrans[,5]),xlim=c(-axis,axis),ylim=c(-axis,axis),col=rgb(t(col2rgb("purple")),alpha=100,maxColorValue=255),pch=19,cex=0.25);

points(log2(cis[,2]/cis[,3]),log2(cis[,4]/cis[,5]),xlim=c(-axis,axis),ylim=c(-axis,axis),col=rgb(t(col2rgb("blue")),alpha=100,maxColorValue=255),pch=19,cex=0.25);

points(log2(trans[,2]/trans[,3]),log2(trans[,4]/trans[,5]),xlim=c(-axis,axis),ylim=c(-axis,axis),col=rgb(t(col2rgb("red")),alpha=100,maxColorValue=255),pch=19,cex=0.25);

dev.off();
















