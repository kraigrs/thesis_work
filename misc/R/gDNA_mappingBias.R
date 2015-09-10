exons <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/zhr_z30.SD_cut20_sepPar2.mosaik.exons.gDNA.txt",header=TRUE,sep="\t");
genes <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/zhr_z30.SD_cut20_sepPar2.mosaik.genes.gDNA.txt",header=TRUE,sep="\t");

hist(log2(exons[,2]/exons[,3]),breaks=50);
hist(log2(genes[,2]/genes[,3]),breaks=50);