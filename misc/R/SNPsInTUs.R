# these paths should be changed, was too lazy to write a script to take args ;)
SNPsPerExon <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/numSNPsPerExon.txt",sep="\t");
SNPsPerGene <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/numSNPsPerGene.txt",sep="\t");

# ~97.5% of exons have 10 or fewer SNPs
sum(SNPsPerExon[,7]<11)/nrow(SNPsPerExon); # = 0.9755268
exon_filtered <- subset(SNPsPerExon,SNPsPerExon[,7]<11);

# ~95% of genes have 100 or fewer SNPs
sum(SNPsPerGene[,5]<101)/nrow(SNPsPerGene); # = 0.9492316
gene_filtered <- subset(SNPsPerGene,SNPsPerGene[,5]<101);

par(mfrow = c(1,2));
hist(exon_filtered[,7],main="~97.5% (53055/54386) of exons with 10 or fewer SNPs",xlab="# SNPs in an exon",breaks=10,freq=FALSE,ylim=c(0,1));
hist(gene_filtered[,5],main="~95% (6918/7288) of genes with 100 or fewer SNPs",xlab="# SNPs in a gene",breaks=100,freq=FALSE,ylim=c(0,1));