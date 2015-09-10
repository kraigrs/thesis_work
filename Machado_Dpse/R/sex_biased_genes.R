############################################################################################
# script to find male-biased genes
############################################################################################

FET_misexpression <- function(locus)
{
	table <- matrix(c(locus[1],locus[3]-locus[1],locus[2],locus[4]-locus[2]),
					nrow=2,
					dimnames=list(c("gene","not gene"),c("sample1","sample2")));
	test <- fisher.test(table,or=1,alternative="two.sided",conf.level=0.95);
	result <- c(test$estimate,test$conf.int[1],test$conf.int[2],test$p.value);
	return(result);
}

create_table <- function(mat)
{
	chrXL <- sum(table(mat$chromosome)[grep("^XL",names(table(mat$chromosome)),perl=TRUE)]);
	chr4 <- sum(table(mat$chromosome)[grep("^4",names(table(mat$chromosome)),perl=TRUE)]);
	chr3 <- sum(table(mat$chromosome)[grep("^3",names(table(mat$chromosome)),perl=TRUE)]);
	chrXR <- sum(table(mat$chromosome)[grep("^XR",names(table(mat$chromosome)),perl=TRUE)]);
	chr2 <- sum(table(mat$chromosome)[grep("^2",names(table(mat$chromosome)),perl=TRUE)]);

    val <- c(chrXL,chr4,chr3,chrXR,chr2);
    return(val);
}

############################################################################################

#F_carcass <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_F_carcass_reg_div.txt",header=TRUE,sep="\t");
#M_carcass <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_M_carcass_reg_div.txt",header=TRUE,sep="\t");

#ovaries <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_ovaries_reg_div.txt",header=TRUE,sep="\t");
#testes <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_testes_reg_div.txt",header=TRUE,sep="\t");

F_carcass <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_F_carcass_MOI.txt",header=TRUE,sep="\t");
M_carcass <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_M_carcass_MOI.txt",header=TRUE,sep="\t");

ovaries <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_ovaries_MOI.txt",header=TRUE,sep="\t");
testes <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_testes_MOI.txt",header=TRUE,sep="\t");

carcass <- merge(M_carcass,F_carcass,by.x="gene",by.y="gene",suffixes=c(".M",".F"));
gonads <- merge(testes,ovaries,by.x="gene",by.y="gene",suffixes=c(".M",".F"));

list <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/genomes/dpse_r3.1.genes.bed",header=FALSE,sep="\t");
colnames(list) <- c("chromosome","start","stop","gene","empty","strand");

FBid2GLEANR <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/genomes/FBid2GLEANR.txt",header=TRUE,sep="\t");

annotation <- merge(FBid2GLEANR,list,by.x="gene",by.y="gene");

# sex-biased genes according to Assis et al.

Assis_MBG <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/Assis_Dpse_MBG.txt",header=TRUE,sep="\t");
Assis_FBG <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/Assis_Dpse_FBG.txt",header=TRUE,sep="\t");

sex_tissue_MBG <- Assis_MBG[which(Assis_MBG$Sex.Tissue.specific == 1),1];
sex_tissue_FBG <- Assis_FBG[which(Assis_FBG$Sex.Tissue.specific == 1),1];

# sex-biased genes according to Zhang et al.

Zhang_MBG <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/Zhang_Dpse_MBG.txt",header=TRUE,sep="\t");
Zhang_FBG <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/Zhang_Dpse_FBG.txt",header=TRUE,sep="\t");

temp <- merge(annotation,Zhang_MBG,by.x="GLEANR",by.y="GleanR.ID");
sex_biased_testes <- merge(testes,temp,by.x="gene",by.y="gene");

temp <- merge(annotation,Zhang_FBG,by.x="GLEANR",by.y="GleanR.ID");
sex_biased_ovaries <- merge(ovaries,temp,by.x="gene",by.y="gene");

###########################
# sex-biased genes by FET #
###########################

sex_bias_carcass_results <- t(apply(as.matrix(cbind(
carcass[,2]+carcass[,3]+carcass[,4],
carcass[,6]+carcass[,7]+carcass[,8],
rep(sum(carcass[,2:4],nrow(carcass))),
rep(sum(carcass[,6:8],nrow(carcass))))),
1,FET_misexpression));

qval <- p.adjust(sex_bias_carcass_results[,4],method="fdr");
OR <- sex_bias_carcass_results[,1];
sex_bias_carcass <- cbind(carcass,OR,qval);

sex_bias_gonads_results <- t(apply(as.matrix(cbind(
gonads[,2]+gonads[,3]+gonads[,4],
gonads[,6]+gonads[,7]+gonads[,8],
rep(sum(gonads[,2:4],nrow(gonads))),
rep(sum(gonads[,6:8],nrow(gonads))))),
1,FET_misexpression));

qval <- p.adjust(sex_bias_gonads_results[,4],method="fdr");
OR <- sex_bias_gonads_results[,1];
sex_bias_gonads <- cbind(gonads,OR,qval);

###################################
# sex-biased genes by fold cutoff #
###################################

fold <- 2;

# carcass

fold_change <- log(((carcass[,2]+carcass[,3]+carcass[,4])/sum(carcass[,2:4]))/((carcass[,6]+carcass[,7]+carcass[,8])/sum(carcass[,6:8])),base=fold);
sex_bias_carcass <- cbind(carcass,fold_change);

MBG <- merge(sex_bias_carcass[which(sex_bias_carcass$fold_change > 1),],list,by.x="gene",by.y="gene");
FBG <- merge(sex_bias_carcass[which(sex_bias_carcass$fold_change < -1),],list,by.x="gene",by.y="gene");
not <- merge(sex_bias_carcass[which(sex_bias_carcass$fold_change >= -1 & sex_bias_carcass$fold_change <= 1),],list,by.x="gene",by.y="gene");

write.table(MBG[,1:10],file="/Users/kraigrs/Wittkopp/Machado_Dpse/fold_change_carcass_MBG.txt",quote=FALSE,sep="\t",row.names=FALSE);
write.table(FBG[,1:10],file="/Users/kraigrs/Wittkopp/Machado_Dpse/fold_change_carcass_FBG.txt",quote=FALSE,sep="\t",row.names=FALSE);

fisher.test(rbind(create_table(MBG),create_table(not)),workspace=100000000);
fisher.test(rbind(create_table(FBG),create_table(not)),workspace=100000000);

# gonads

fold_change <- log(((gonads[,2]+gonads[,3]+gonads[,4])/sum(gonads[,2:4]))/((gonads[,6]+gonads[,7]+gonads[,8])/sum(gonads[,6:8])),base=fold);
sex_bias_gonads <- cbind(gonads,fold_change);

MBG <- merge(sex_bias_gonads[which(sex_bias_gonads$fold_change > 1),],list,by.x="gene",by.y="gene");
FBG <- merge(sex_bias_gonads[which(sex_bias_gonads$fold_change < -1),],list,by.x="gene",by.y="gene");
not <- merge(sex_bias_gonads[which(sex_bias_gonads$fold_change >= -1 & sex_bias_gonads$fold_change <= 1),],list,by.x="gene",by.y="gene");

write.table(MBG[,1:10],file="/Users/kraigrs/Wittkopp/Machado_Dpse/fold_change_gonads_MBG.txt",quote=FALSE,sep="\t",row.names=FALSE);
write.table(FBG[,1:10],file="/Users/kraigrs/Wittkopp/Machado_Dpse/fold_change_gonads_FBG.txt",quote=FALSE,sep="\t",row.names=FALSE);

fisher.test(rbind(create_table(MBG),create_table(not)),workspace=100000000);
fisher.test(rbind(create_table(FBG),create_table(not)),workspace=100000000);

##################################
# find sex tissue-specific genes #
##################################

i <- testes$gene %in% gonads$gene;
MBG <- merge(testes[!i,],list,by.x="gene",by.y="gene");

j <- ovaries$gene %in% gonads$gene;
FBG <- merge(ovaries[!j,],list,by.x="gene",by.y="gene");

k <- MBG$gene %in% M_carcass$gene;
testis_specific <- MBG[!k,];

l <- FBG$gene %in% F_carcass$gene;
ovary_specific <- FBG[!l,];

write.table(testis_specific[,1:5],file="/Users/kraigrs/Wittkopp/Machado_Dpse/testis_specific.txt",quote=FALSE,sep="\t",row.names=FALSE);

write.table(ovary_specific[,1:5],file="/Users/kraigrs/Wittkopp/Machado_Dpse/ovary_specific.txt",quote=FALSE,sep="\t",row.names=FALSE);

plot(log2(gonads$par1.testes + gonads$par2.testes + gonads$hyb.testes),
     log2(gonads$par1.ovaries + gonads$par2.ovaries + gonads$hyb.ovaries),
     xlim=c(log2(60),log2(400000)),ylim=c(log2(60),log2(400000)),
     xlab="log2(testes)",ylab="log2(ovaries)",
     pch=16,col=rgb(0,0,0,0.5),cex=0.5);
abline(a=0,b=1,col="red");

plot((log2(gonads$par1.testes+gonads$par2.testes+gonads$hyb.testes)+log2(gonads$par1.ovaries+gonads$par2.ovaries+gonads$hyb.ovaries))/2,
     log2((gonads$par1.testes+gonads$par2.testes+gonads$hyb.testes)/(gonads$par1.ovaries+gonads$par2.ovaries+gonads$hyb.ovaries)),
     ylim=c(-10,10),xlab="avg. log2(expression)",ylab="log2(testes/ovaries)",
     pch=16,col=rgb(0,0,0,0.5),cex=0.5);


