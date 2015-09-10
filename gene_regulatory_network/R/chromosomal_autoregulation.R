########## expression data ##########

chrom_auto <- read.table("/Users/kraigrs/Wittkopp/gene_regulatory_network/Dmel_chromosomal_autoregulation.txt",header=TRUE,sep="\t");

chrom_auto <- chrom_auto[which(chrom_auto$chromosome != 4),];

chrom_auto$chromosome <- factor(chrom_auto$chromosome,levels=c("X","2L","2R","3L","3R"));

boxplot(autoreg/targets ~ chromosome,data=chrom_auto,outline=FALSE,varwidth=TRUE,ylim=c(0,1),xlab="chromosome",ylab="proportion autoregulated targets",main="regulators");

boxplot(targets ~ chromosome,data=chrom_auto,outline=FALSE,varwidth=TRUE,xlab="chromosome",ylab="# targets",main="regulators");

########## Kellis data ##########

chrom_auto <- read.table("/Users/kraigrs/Wittkopp/gene_regulatory_network/Kellis_chromosomal_autoregulation.txt",header=TRUE,sep="\t");

chrom_auto <- chrom_auto[which(chrom_auto$chromosome != 4),];

chrom_auto$chromosome <- factor(chrom_auto$chromosome,levels=c("X","2L","2R","3L","3R"));

boxplot(autoreg/targets ~ chromosome,data=chrom_auto,outline=FALSE,varwidth=TRUE,ylim=c(0,1),xlab="chromosome",ylab="proportion autoregulated targets",main="regulators (Kellis)");

boxplot(targets ~ chromosome,data=chrom_auto,outline=FALSE,varwidth=TRUE,xlab="chromosome",ylab="# targets targets",main="regulators");

########## human data (raw) ##########

chrom_auto <- read.table("/Users/kraigrs/Wittkopp/gene_regulatory_network/enets1.Proximal_raw.autoreg.txt",header=TRUE,sep="\t");

chrom_auto <- chrom_auto[grep("^[0-9X]",chrom_auto$chromosome,perl=TRUE),];

chrom_auto$chromosome <- factor(chrom_auto$chromosome,levels=c("X",seq(1,22,1)));

boxplot(autoreg/targets ~ chromosome,data=chrom_auto,varwidth=TRUE,ylim=c(0,0.2),xlab="chromosome",ylab="proportion targets on same chromosome targets",main="regulators (human data)");

boxplot(targets ~ chromosome,data=chrom_auto,varwidth=TRUE,xlab="chromosome",ylab="# targets",main="regulators");

########## human data (filtered) ##########

chrom_auto <- read.table("/Users/kraigrs/Wittkopp/gene_regulatory_network/enets2.Proximal_filtered.autoreg.txt",header=TRUE,sep="\t");

chrom_auto <- chrom_auto[grep("^[0-9X]",chrom_auto$chromosome,perl=TRUE),];

chrom_auto$chromosome <- factor(chrom_auto$chromosome,levels=c("X",seq(1,22,1)));

boxplot(autoreg/targets ~ chromosome,data=chrom_auto,varwidth=TRUE,ylim=c(0,0.4),xlab="chromosome",ylab="proportion targets on same chromosome targets",main="regulators (human data)");

boxplot(targets ~ chromosome,data=chrom_auto,varwidth=TRUE,xlab="chromosome",ylab="# targets",main="regulators");