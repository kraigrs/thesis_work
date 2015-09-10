data <- read.table("/Users/kraigrs/Wittkopp/gene_regulatory_network/regulator_summary.txt",header=TRUE,sep="\t");

data$reg_div <- factor(data$reg_div,levels=c("cis","trans","cis+trans","cisXtrans","compensatory","conserved","ambiguous"));

boxplot(trans/FB_targets ~ reg_div, data=data, varwidth=TRUE,ylab="proportion of targets with trans-acting effects");