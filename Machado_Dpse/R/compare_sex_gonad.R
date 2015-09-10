############################################################################################
# script to compare sexes
############################################################################################

# carcass

data1 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_F_carcass_reg_div.txt",header=TRUE,sep="\t");
data2 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_M_carcass_reg_div_no_X.txt",header=TRUE,sep="\t");

data <- merge(data1,data2,by.x="gene",by.y="gene");
final <- subset(data,!is.na(percent_cis.x) & !is.na(percent_cis.y));

wilcox.test(final$percent_cis.x*100,final$percent_cis.y*100,paired=FALSE,alternative="two.sided",conf.int=TRUE);

pdf(file="/Users/kraigrs/Wittkopp/Machado_Dpse/plots/carcass_percent_cis.pdf");
boxplot(final$percent_cis.x*100,at=1,ylim=c(0,100),xlim=c(0.5,2));
boxplot(final$percent_cis.y*100,add=TRUE,at=1.5,ylim=c(0,100),xlim=c(0.5,2),ylab="% cis",main="carcass");
dev.off();


# gonads

data1 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_ovaries_reg_div.txt",header=TRUE,sep="\t");
data2 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_testes_reg_div_no_X.txt",header=TRUE,sep="\t");

data <- merge(data1,data2,by.x="gene",by.y="gene");
final <- subset(data,!is.na(percent_cis.x) & !is.na(percent_cis.y));

wilcox.test(final$percent_cis.x*100,final$percent_cis.y*100,paired=FALSE,alternative="two.sided",conf.int=TRUE);

pdf(file="/Users/kraigrs/Wittkopp/Machado_Dpse/plots/gonad_percent_cis.pdf");
boxplot(final$percent_cis.x*100,at=1,ylim=c(0,100),xlim=c(0.5,2));
boxplot(final$percent_cis.y*100,add=TRUE,at=1.5,ylim=c(0,100),xlim=c(0.5,2),ylab="% cis",main="gonad");
dev.off();
