############################################################################################
# script to calculate rho values among sex tissue-specific and sex-biased genes
############################################################################################

# functions

boot <- function(mat)
{
	set1 <- sample(seq(1,nrow(mat),1),nrow(mat),replace=TRUE);
	set2 <- sample(seq(1,nrow(mat),1),nrow(mat),replace=TRUE);
	val1 <- 1-cor(log10(mat$Dpse_TL[set1]),log10(mat$Dbog_Toro1[set1]),use="pairwise.complete.obs",method="spearman");
	val2 <- 1-cor(log10(mat$Hyb_Dpse[set2]),log10(mat$Hyb_Dbog[set2]),use="pairwise.complete.obs",method="spearman");
	return(c(val1,val2));
}

############################################################################################

# allele-specific parental and hybrid expression

# H6
#set.seed(12345); # F_carcass
#set.seed(67891); # ovaries
#set.seed(34567); # M_carcass
set.seed(89123); # testes

# H5
#set.seed(23456); # F_carcass
#set.seed(78912); # ovaries
#set.seed(45678); # M_carcass
#set.seed(91234); # testes

#data <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_F_carcass_reg_div.txt",header=TRUE,sep="\t");
#data <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_ovaries_reg_div.txt",header=TRUE,sep="\t");
#data <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_M_carcass_reg_div.txt",header=TRUE,sep="\t");
data <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_testes_reg_div.txt",header=TRUE,sep="\t");

list <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/genomes/dpse_r3.1.genes.bed",header=FALSE,sep="\t");
colnames(list) <- c("chromosome","start","stop","gene","empty","strand");

data_annotation <- merge(data,list,by.x="gene",by.y="gene");

##### sex-biased genes #####

#MBG <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/fold_change_carcass_MBG.txt",header=TRUE,sep="\t");
#FBG <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/fold_change_carcass_FBG.txt",header=TRUE,sep="\t");

#MBG <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/fold_change_gonads_MBG.txt",header=TRUE,sep="\t");
#FBG <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/fold_change_gonads_FBG.txt",header=TRUE,sep="\t");

MBG <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/testis_specific.txt",header=TRUE,sep="\t");
FBG <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/ovary_specific.txt",header=TRUE,sep="\t");

male <- merge(data_annotation,MBG,by.x="gene",by.y="gene");

male_XL <- male[grep("^XL",male$chromosome,perl=TRUE),];
male_4 <- male[grep("^4",male$chromosome,perl=TRUE),];
male_3 <- male[grep("^3",male$chromosome,perl=TRUE),];
male_XR <- male[grep("^XR",male$chromosome,perl=TRUE),];
male_2 <- male[grep("^2",male$chromosome,perl=TRUE),];
male_all <- male;
male_auto <- rbind(male_4,male_3,male_2);
male_X <- rbind(male_XL,male_XR);

female <- merge(data_annotation,FBG,by.x="gene",by.y="gene");

female_XL <- female[grep("^XL",female$chromosome,perl=TRUE),];
female_4 <- female[grep("^4",female$chromosome,perl=TRUE),];
female_3 <- female[grep("^3",female$chromosome,perl=TRUE),];
female_XR <- female[grep("^XR",female$chromosome,perl=TRUE),];
female_2 <- female[grep("^2",female$chromosome,perl=TRUE),];
female_all <- female;
female_auto <- rbind(female_4,female_3,female_2);
female_X <- rbind(female_XL,female_XR);

genes <- merge(male,female,all=TRUE)[,1];

index <- data_annotation$gene %in% genes;
data <- data_annotation[!index,];

data_XL <- data[grep("^XL",data$chromosome,perl=TRUE),];
data_4 <- data[grep("^4",data$chromosome,perl=TRUE),];
data_3 <- data[grep("^3",data$chromosome,perl=TRUE),];
data_XR <- data[grep("^XR",data$chromosome,perl=TRUE),];
data_2 <- data[grep("^2",data$chromosome,perl=TRUE),];
data_all <- data;
data_auto <- rbind(data_4,data_3,data_2);
data_X <- rbind(data_XL,data_XR);

# bootstrap

male_XL_boot <- t(replicate(10000,boot(male_XL)));
#male_4_boot <- t(replicate(10000,boot(male_4)));
#male_3_boot <- t(replicate(10000,boot(male_3)));
male_XR_boot <- t(replicate(10000,boot(male_XR)));
#male_2_boot <- t(replicate(10000,boot(male_2)));
male_all_boot <- t(replicate(10000,boot(male_all)));
male_auto_boot <- t(replicate(10000,boot(male_auto)));
male_X_boot <- t(replicate(10000,boot(male_X)));

female_XL_boot <- t(replicate(10000,boot(female_XL)));
#female_4_boot <- t(replicate(10000,boot(female_4)));
#female_3_boot <- t(replicate(10000,boot(female_3)));
female_XR_boot <- t(replicate(10000,boot(female_XR)));
#female_2_boot <- t(replicate(10000,boot(female_2)));
female_all_boot <- t(replicate(10000,boot(female_all)));
female_auto_boot <- t(replicate(10000,boot(female_auto)));
female_X_boot <- t(replicate(10000,boot(female_X)));

chr_XL_boot <- t(replicate(10000,boot(data_XL)));
#chr_4_boot <- t(replicate(10000,boot(data_4)));
#chr_3_boot <- t(replicate(10000,boot(data_3)));
chr_XR_boot <- t(replicate(10000,boot(data_XR)));
#chr_2_boot <- t(replicate(10000,boot(data_2)));
chr_all_boot <- t(replicate(10000,boot(data_all)));
chr_auto_boot <- t(replicate(10000,boot(data_auto)));
chr_X_boot <- t(replicate(10000,boot(data_X)));

# parent1 to parent2

#pdf(file="/Users/kraigrs/Wittkopp/Machado_Dpse/plots/TL_Toro1_H6_F_carcass_tot_expr_div_sex_bias.pdf",height=5,width=5);
#pdf(file="/Users/kraigrs/Wittkopp/Machado_Dpse/plots/TL_Toro1_H6_ovaries_tot_expr_div_sex_bias.pdf",height=5,width=5);
#pdf(file="/Users/kraigrs/Wittkopp/Machado_Dpse/plots/TL_Toro1_H6_M_carcass_tot_expr_div_sex_bias.pdf",height=5,width=5);
#pdf(file="/Users/kraigrs/Wittkopp/Machado_Dpse/plots/TL_Toro1_H6_testes_tot_expr_div_sex_bias.pdf",height=5,width=5);

#pdf(file="/Users/kraigrs/Wittkopp/Machado_Dpse/plots/TL_Toro1_H6_ovaries_tot_expr_div_ovary_specific.pdf",height=5,width=5);
pdf(file="/Users/kraigrs/Wittkopp/Machado_Dpse/plots/TL_Toro1_H6_testes_tot_expr_div_testis_specific.pdf",height=5,width=5);

maleXL <- 1-cor(log10(male_XL$Dpse_TL),log10(male_XL$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");
#male4 <- 1-cor(log10(male_4$Dpse_TL),log10(male_4$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");
#male3 <- 1-cor(log10(male_3$Dpse_TL),log10(male_3$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");
maleXR <- 1-cor(log10(male_XR$Dpse_TL),log10(male_XR$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");
#male2 <- 1-cor(log10(male_2$Dpse_TL),log10(male_2$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");
male_overall <- 1-cor(log10(male_all$Dpse_TL),log10(male_all$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");
m_auto <- 1-cor(log10(male_auto$Dpse_TL),log10(male_auto$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");
maleX <- 1-cor(log10(male_X$Dpse_TL),log10(male_X$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");

femaleXL <- 1-cor(log10(female_XL$Dpse_TL),log10(female_XL$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");
#female4 <- 1-cor(log10(female_4$Dpse_TL),log10(female_4$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");
#female3 <- 1-cor(log10(female_3$Dpse_TL),log10(female_3$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");
femaleXR <- 1-cor(log10(female_XR$Dpse_TL),log10(female_XR$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");
#female2 <- 1-cor(log10(female_2$Dpse_TL),log10(female_2$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");
female_overall <- 1-cor(log10(female_all$Dpse_TL),log10(female_all$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");
f_auto <- 1-cor(log10(female_auto$Dpse_TL),log10(female_auto$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");
femaleX <- 1-cor(log10(female_X$Dpse_TL),log10(female_X$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");

chrXL <- 1-cor(log10(data_XL$Dpse_TL),log10(data_XL$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");
#chr4 <- 1-cor(log10(data_4$Dpse_TL),log10(data_4$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");
#chr3 <- 1-cor(log10(data_3$Dpse_TL),log10(data_3$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");
chrXR <- 1-cor(log10(data_XR$Dpse_TL),log10(data_XR$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");
#chr2 <- 1-cor(log10(data_2$Dpse_TL),log10(data_2$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");
overall <- 1-cor(log10(data_all$Dpse_TL),log10(data_all$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");
auto <- 1-cor(log10(data_auto$Dpse_TL),log10(data_auto$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");
chrX <- 1-cor(log10(data_X$Dpse_TL),log10(data_X$Dbog_Toro1),use="pairwise.complete.obs",method="spearman");

plot(seq(0,20,1),rep(0,21),type="n",ylim=c(0,0.4),ylab="1-rho",xlab="",xaxt="n",main="parent1 vs. parent2");

# X

points(1,maleX,pch=16,col="cornflowerblue"); points(c(1,1),quantile(male_X_boot[,1],probs=c(0.025,0.975)),pch="-",col="cornflowerblue");
lines(c(1,1),quantile(male_X_boot[,1],probs=c(0.025,0.975)),lty=1,col="cornflowerblue");
points(2,chrX,pch=16,col="gray"); points(c(2,2),quantile(chr_X_boot[,1],probs=c(0.025,0.975)),pch="-",col="gray");
lines(c(2,2),quantile(chr_X_boot[,1],probs=c(0.025,0.975)),lty=1,col="gray");
points(3,femaleX,pch=16,col="plum1"); points(c(3,3),quantile(female_X_boot[,1],probs=c(0.025,0.975)),pch="-",col="plum1");
lines(c(3,3),quantile(female_X_boot[,1],probs=c(0.025,0.975)),lty=1,col="plum1");

# auto

points(5,m_auto,pch=16,col="cornflowerblue"); points(c(5,5),quantile(male_auto_boot[,1],probs=c(0.025,0.975)),pch="-",col="cornflowerblue");
lines(c(5,5),quantile(male_auto_boot[,1],probs=c(0.025,0.975)),lty=1,col="cornflowerblue");
points(6,auto,pch=16,col="gray"); points(c(6,6),quantile(chr_auto_boot[,1],probs=c(0.025,0.975)),pch="-",col="gray");
lines(c(6,6),quantile(chr_auto_boot[,1],probs=c(0.025,0.975)),lty=1,col="gray");
points(7,f_auto,pch=16,col="plum1"); points(c(7,7),quantile(female_auto_boot[,1],probs=c(0.025,0.975)),pch="-",col="plum1");
lines(c(7,7),quantile(female_auto_boot[,1],probs=c(0.025,0.975)),lty=1,col="plum1");

# XL

points(9,maleXL,pch=16,col="cornflowerblue"); points(c(9,9),quantile(male_XL_boot[,1],probs=c(0.025,0.975)),pch="-",col="cornflowerblue");
lines(c(9,9),quantile(male_XL_boot[,1],probs=c(0.025,0.975)),lty=1,col="cornflowerblue");
points(10,chrXL,pch=16,col="gray"); points(c(10,10),quantile(chr_XL_boot[,1],probs=c(0.025,0.975)),pch="-",col="gray");
lines(c(10,10),quantile(chr_XL_boot[,1],probs=c(0.025,0.975)),lty=1,col="gray");
points(11,femaleXL,pch=16,col="plum1"); points(c(11,11),quantile(female_XL_boot[,1],probs=c(0.025,0.975)),pch="-",col="plum1");
lines(c(11,11),quantile(female_XL_boot[,1],probs=c(0.025,0.975)),lty=1,col="plum1");

# XR

points(13,maleXR,pch=16,col="cornflowerblue"); points(c(13,13),quantile(male_XR_boot[,1],probs=c(0.025,0.975)),pch="-",col="cornflowerblue");
lines(c(13,13),quantile(male_XR_boot[,1],probs=c(0.025,0.975)),lty=1,col="cornflowerblue");
points(14,chrXR,pch=16,col="gray"); points(c(14,14),quantile(chr_XR_boot[,1],probs=c(0.025,0.975)),pch="-",col="gray");
lines(c(14,14),quantile(chr_XR_boot[,1],probs=c(0.025,0.975)),lty=1,col="gray");
points(15,femaleXR,pch=16,col="plum1"); points(c(15,15),quantile(female_XR_boot[,1],probs=c(0.025,0.975)),pch="-",col="plum1");
lines(c(15,15),quantile(female_XR_boot[,1],probs=c(0.025,0.975)),lty=1,col="plum1");

# all

points(17,male_overall,pch=16,col="cornflowerblue"); points(c(17,17),quantile(male_all_boot[,1],probs=c(0.025,0.975)),pch="-",col="cornflowerblue");
lines(c(17,17),quantile(male_all_boot[,1],probs=c(0.025,0.975)),lty=1,col="cornflowerblue");
points(18,overall,pch=16,col="gray"); points(c(18,18),quantile(chr_all_boot[,1],probs=c(0.025,0.975)),pch="-",col="gray");
lines(c(18,18),quantile(chr_all_boot[,1],probs=c(0.025,0.975)),lty=1,col="gray");
points(19,female_overall,pch=16,col="plum1"); points(c(19,19),quantile(female_all_boot[,1],probs=c(0.025,0.975)),pch="-",col="plum1");
lines(c(19,19),quantile(female_all_boot[,1],probs=c(0.025,0.975)),lty=1,col="plum1");

axis(1,at=c(2,6,10,14,18),labels=c("XLR","auto","X","neo-X","all"));

abline(v=c(8,16),col="black");

dev.off();

# ASE in hybrid

#pdf(file="/Users/kraigrs/Wittkopp/Machado_Dpse/plots/TL_Toro1_H6_F_carcass_cis_expr_div_sex_bias.pdf",height=5,width=5);
#pdf(file="/Users/kraigrs/Wittkopp/Machado_Dpse/plots/TL_Toro1_H6_ovaries_cis_expr_div_sex_bias.pdf",height=5,width=5);

pdf(file="/Users/kraigrs/Wittkopp/Machado_Dpse/plots/TL_Toro1_H6_ovaries_cis_expr_div_ovary_specific.pdf",height=5,width=5);


maleXL <- 1-cor(log10(male_XL$Hyb_Dpse),log10(male_XL$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");
#male4 <- 1-cor(log10(male_4$Hyb_Dpse),log10(male_4$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");
#male3 <- 1-cor(log10(male_3$Hyb_Dpse),log10(male_3$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");
maleXR <- 1-cor(log10(male_XR$Hyb_Dpse),log10(male_XR$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");
#male2 <- 1-cor(log10(male_2$Hyb_Dpse),log10(male_2$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");
male_overall <- 1-cor(log10(male_all$Hyb_Dpse),log10(male_all$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");
m_auto <- 1-cor(log10(male_auto$Hyb_Dpse),log10(male_auto$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");
maleX <- 1-cor(log10(male_X$Hyb_Dpse),log10(male_X$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");

femaleXL <- 1-cor(log10(female_XL$Hyb_Dpse),log10(female_XL$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");
#female4 <- 1-cor(log10(female_4$Hyb_Dpse),log10(female_4$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");
#female3 <- 1-cor(log10(female_3$Hyb_Dpse),log10(female_3$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");
femaleXR <- 1-cor(log10(female_XR$Hyb_Dpse),log10(female_XR$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");
#female2 <- 1-cor(log10(female_2$Hyb_Dpse),log10(female_2$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");
female_overall <- 1-cor(log10(female_all$Hyb_Dpse),log10(female_all$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");
f_auto <- 1-cor(log10(female_auto$Hyb_Dpse),log10(female_auto$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");
femaleX <- 1-cor(log10(female_X$Hyb_Dpse),log10(female_X$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");

chrXL <- 1-cor(log10(data_XL$Hyb_Dpse),log10(data_XL$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");
#chr4 <- 1-cor(log10(data_4$Hyb_Dpse),log10(data_4$Hyb_Dbog)d,use="pairwise.complete.obs",method="spearman");
#chr3 <- 1-cor(log10(data_3$Hyb_Dpse),log10(data_3$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");
chrXR <- 1-cor(log10(data_XR$Hyb_Dpse),log10(data_XR$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");
#chr2 <- 1-cor(log10(data_2$Hyb_Dpse),log10(data_2$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");
overall <- 1-cor(log10(data_all$Hyb_Dpse),log10(data_all$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");
auto <- 1-cor(log10(data_auto$Hyb_Dpse),log10(data_auto$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");
chrX <- 1-cor(log10(data_X$Hyb_Dpse),log10(data_X$Hyb_Dbog),use="pairwise.complete.obs",method="spearman");

plot(seq(0,20,1),rep(0,21),type="n",ylim=c(0,0.4),ylab="1-rho",xlab="",xaxt="n",main="hybrid allele1 vs. allele2");

# X

points(1,maleX,pch=16,col="cornflowerblue"); points(c(1,1),quantile(male_X_boot[,2],probs=c(0.025,0.975)),pch="-",col="cornflowerblue");
lines(c(1,1),quantile(male_X_boot[,2],probs=c(0.025,0.975)),lty=1,col="cornflowerblue");
points(2,chrX,pch=16,col="gray"); points(c(2,2),quantile(chr_X_boot[,2],probs=c(0.025,0.975)),pch="-",col="gray");
lines(c(2,2),quantile(chr_X_boot[,2],probs=c(0.025,0.975)),lty=1,col="gray");
points(3,femaleX,pch=16,col="plum1"); points(c(3,3),quantile(female_X_boot[,2],probs=c(0.025,0.975)),pch="-",col="plum1");
lines(c(3,3),quantile(female_X_boot[,2],probs=c(0.025,0.975)),lty=1,col="plum1");

# auto

points(5,m_auto,pch=16,col="cornflowerblue"); points(c(5,5),quantile(male_auto_boot[,2],probs=c(0.025,0.975)),pch="-",col="cornflowerblue");
lines(c(5,5),quantile(male_auto_boot[,2],probs=c(0.025,0.975)),lty=1,col="cornflowerblue");
points(6,auto,pch=16,col="gray"); points(c(6,6),quantile(chr_auto_boot[,2],probs=c(0.025,0.975)),pch="-",col="gray");
lines(c(6,6),quantile(chr_auto_boot[,2],probs=c(0.025,0.975)),lty=1,col="gray");
points(7,f_auto,pch=16,col="plum1"); points(c(7,7),quantile(female_auto_boot[,2],probs=c(0.025,0.975)),pch="-",col="plum1");
lines(c(7,7),quantile(female_auto_boot[,2],probs=c(0.025,0.975)),lty=1,col="plum1");

# XL

points(9,maleXL,pch=16,col="cornflowerblue"); points(c(9,9),quantile(male_XL_boot[,2],probs=c(0.025,0.975)),pch="-",col="cornflowerblue");
lines(c(9,9),quantile(male_XL_boot[,2],probs=c(0.025,0.975)),lty=1,col="cornflowerblue");
points(10,chrXL,pch=16,col="gray"); points(c(10,10),quantile(chr_XL_boot[,2],probs=c(0.025,0.975)),pch="-",col="gray");
lines(c(10,10),quantile(chr_XL_boot[,2],probs=c(0.025,0.975)),lty=1,col="gray");
points(11,femaleXL,pch=16,col="plum1"); points(c(11,11),quantile(female_XL_boot[,2],probs=c(0.025,0.975)),pch="-",col="plum1");
lines(c(11,11),quantile(female_XL_boot[,2],probs=c(0.025,0.975)),lty=1,col="plum1");

# XR

points(13,maleXR,pch=16,col="cornflowerblue"); points(c(13,13),quantile(male_XR_boot[,2],probs=c(0.025,0.975)),pch="-",col="cornflowerblue");
lines(c(13,13),quantile(male_XR_boot[,2],probs=c(0.025,0.975)),lty=1,col="cornflowerblue");
points(14,chrXR,pch=16,col="gray"); points(c(14,14),quantile(chr_XR_boot[,2],probs=c(0.025,0.975)),pch="-",col="gray");
lines(c(14,14),quantile(chr_XR_boot[,2],probs=c(0.025,0.975)),lty=1,col="gray");
points(15,femaleXR,pch=16,col="plum1"); points(c(15,15),quantile(female_XR_boot[,2],probs=c(0.025,0.975)),pch="-",col="plum1");
lines(c(15,15),quantile(female_XR_boot[,2],probs=c(0.025,0.975)),lty=1,col="plum1");

# all

points(17,male_overall,pch=16,col="cornflowerblue"); points(c(17,17),quantile(male_all_boot[,2],probs=c(0.025,0.975)),pch="-",col="cornflowerblue");
lines(c(17,17),quantile(male_all_boot[,2],probs=c(0.025,0.975)),lty=1,col="cornflowerblue");
points(18,overall,pch=16,col="gray"); points(c(18,18),quantile(chr_all_boot[,2],probs=c(0.025,0.975)),pch="-",col="gray");
lines(c(18,18),quantile(chr_all_boot[,2],probs=c(0.025,0.975)),lty=1,col="gray");
points(19,female_overall,pch=16,col="plum1"); points(c(19,19),quantile(female_all_boot[,2],probs=c(0.025,0.975)),pch="-",col="plum1");
lines(c(19,19),quantile(female_all_boot[,2],probs=c(0.025,0.975)),lty=1,col="plum1");

axis(1,at=c(2,6,10,14,18),labels=c("XLR","auto","X","neo-X","all"));

abline(v=c(8,16),col="black");

dev.off();



####################################################################
# old stuff

plot(seq(0,20,1),rep(male_overall,21),type="l",lty=1,ylim=c(0,0.4),ylab="1-rho",xlab="",xaxt="n",main="parent1 vs. parent2",col="cornflowerblue");
lines(c(0,20),rep(quantile(male_all_boot[,1],probs=0.025),2),col="cornflowerblue",lty=2);
lines(c(0,20),rep(quantile(male_all_boot[,1],probs=0.975),2),col="cornflowerblue",lty=2);

lines(c(0,20),rep(overall,2),col="gray",lty=1);
lines(c(0,20),rep(quantile(chr_all_boot[,1],probs=0.025),2),col="gray",lty=2);
lines(c(0,20),rep(quantile(chr_all_boot[,1],probs=0.975),2),col="gray",lty=2);

lines(c(0,20),rep(female_overall,2),col="plum1",lty=1);
lines(c(0,20),rep(quantile(female_all_boot[,1],probs=0.025),2),col="plum1",lty=2);
lines(c(0,20),rep(quantile(female_all_boot[,1],probs=0.975),2),col="plum1",lty=2);

# chrXL
points(1,maleXL,pch=16,col="cornflowerblue"); points(c(1,1),quantile(male_XL_boot[,1],probs=c(0.025,0.975)),pch="-",col="cornflowerblue");
lines(c(1,1),quantile(male_XL_boot[,1],probs=c(0.025,0.975)),lty=1,col="cornflowerblue");
points(2,chrXL,pch=16,col="gray"); points(c(2,2),quantile(chr_XL_boot[,1],probs=c(0.025,0.975)),pch="-",col="gray");
lines(c(2,2),quantile(chr_XL_boot[,1],probs=c(0.025,0.975)),lty=1,col="gray");
points(3,femaleXL,pch=16,col="plum1"); points(c(3,3),quantile(female_XL_boot[,1],probs=c(0.025,0.975)),pch="-",col="plum1");
lines(c(3,3),quantile(female_XL_boot[,1],probs=c(0.025,0.975)),lty=1,col="plum1");

# chr4
points(5,male4,pch=16,col="cornflowerblue"); points(c(5,5),quantile(male_4_boot[,1],probs=c(0.025,0.975)),pch="-",col="cornflowerblue");
lines(c(5,5),quantile(male_4_boot[,1],probs=c(0.025,0.975)),lty=1,col="cornflowerblue");
points(6,chr4,pch=16,col="gray"); points(c(6,6),quantile(chr_4_boot[,1],probs=c(0.025,0.975)),pch="-",col="gray");
lines(c(6,6),quantile(chr_4_boot[,1],probs=c(0.025,0.975)),lty=1,col="gray");
points(7,female4,pch=16,col="plum1"); points(c(7,7),quantile(female_4_boot[,1],probs=c(0.025,0.975)),pch="-",col="plum1");
lines(c(7,7),quantile(female_4_boot[,1],probs=c(0.025,0.975)),lty=1,col="plum1");

# chr3

points(9,male3,pch=16,col="cornflowerblue"); points(c(9,9),quantile(male_3_boot[,1],probs=c(0.025,0.975)),pch="-",col="cornflowerblue");
lines(c(9,9),quantile(male_3_boot[,1],probs=c(0.025,0.975)),lty=1,col="cornflowerblue");
points(10,chr3,pch=16,col="gray"); points(c(10,10),quantile(chr_3_boot[,1],probs=c(0.025,0.975)),pch="-",col="gray");
lines(c(10,10),quantile(chr_3_boot[,1],probs=c(0.025,0.975)),lty=1,col="gray");
points(11,female3,pch=16,col="plum1"); points(c(11,11),quantile(female_3_boot[,1],probs=c(0.025,0.975)),pch="-",col="plum1");
lines(c(11,11),quantile(female_3_boot[,1],probs=c(0.025,0.975)),lty=1,col="plum1");

# chrXR

points(13,maleXR,pch=16,col="cornflowerblue"); points(c(13,13),quantile(male_XR_boot[,1],probs=c(0.025,0.975)),pch="-",col="cornflowerblue");
lines(c(13,13),quantile(male_XR_boot[,1],probs=c(0.025,0.975)),lty=1,col="cornflowerblue");
points(14,chrXR,pch=16,col="gray"); points(c(14,14),quantile(chr_XR_boot[,1],probs=c(0.025,0.975)),pch="-",col="gray");
lines(c(14,14),quantile(chr_XR_boot[,1],probs=c(0.025,0.975)),lty=1,col="gray");
points(15,femaleXR,pch=16,col="plum1"); points(c(15,15),quantile(female_XR_boot[,1],probs=c(0.025,0.975)),pch="-",col="plum1");
lines(c(15,15),quantile(female_XR_boot[,1],probs=c(0.025,0.975)),lty=1,col="plum1");

# chr2

points(17,male2,pch=16,col="cornflowerblue"); points(c(17,17),quantile(male_2_boot[,1],probs=c(0.025,0.975)),pch="-",col="cornflowerblue");
lines(c(17,17),quantile(male_2_boot[,1],probs=c(0.025,0.975)),lty=1,col="cornflowerblue");
points(18,chr2,pch=16,col="gray"); points(c(18,18),quantile(chr_2_boot[,1],probs=c(0.025,0.975)),pch="-",col="gray");
lines(c(18,18),quantile(chr_2_boot[,1],probs=c(0.025,0.975)),lty=1,col="gray");
points(19,female2,pch=16,col="plum1"); points(c(19,19),quantile(female_2_boot[,1],probs=c(0.025,0.975)),pch="-",col="plum1");
lines(c(19,19),quantile(female_2_boot[,1],probs=c(0.025,0.975)),lty=1,col="plum1");

axis(1,at=c(2,6,10,14,18),labels=c("X","4","3","neo-X","2"));