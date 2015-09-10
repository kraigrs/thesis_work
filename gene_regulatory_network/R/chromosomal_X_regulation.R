########## expression data ##########

chrom_X <- read.table("/Users/kraigrs/Wittkopp/gene_regulatory_network/Dmel_chromosomal_X_regulation.txt",header=TRUE,sep="\t");

chrom_X <- chrom_X[which(chrom_X$chromosome != 4),];

chrom_X$chromosome <- factor(chrom_X$chromosome,levels=c("X","2L","2R","3L","3R"));

par(mfrow=c(1,2));

boxplot(Xreg/targets ~ chromosome,data=chrom_X,varwidth=TRUE,ylim=c(0,1),xlab="chromosome",ylab="proportion targets on X",main="expression data");

#boxplot(targets ~ chromosome,data=chrom_X,outline=FALSE,varwidth=TRUE,xlab="chromosome",ylab="# targets",main="expression data");

########## Kellis data ##########

chrom_X <- read.table("/Users/kraigrs/Wittkopp/gene_regulatory_network/Kellis_chromosomal_X_regulation.txt",header=TRUE,sep="\t");

chrom_X <- chrom_X[which(chrom_X$chromosome != 4),];

chrom_X$chromosome <- factor(chrom_X$chromosome,levels=c("X","2L","2R","3L","3R"));

boxplot(Xreg/targets ~ chromosome,data=chrom_X,varwidth=TRUE,ylim=c(0,1),xlab="chromosome",ylab="proportion targets on X",main="Kellis data");

#boxplot(targets ~ chromosome,data=chrom_X,outline=FALSE,varwidth=TRUE,xlab="chromosome",ylab="# targets targets",main="regulators");