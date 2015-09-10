zhr_exons <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/Resequencing/zhr/zhr.gDNA.mosaik.exonout.txt",header=TRUE);
zhr_genes <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/Resequencing/zhr/zhr.gDNA.mosaik.geneout.txt",header=TRUE);

z30_exons <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/Resequencing/z30/z30.gDNA.mosaik.exonout.txt",header=TRUE);
z30_genes <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/Resequencing/z30/z30.gDNA.mosaik.geneout.txt",header=TRUE);

obj <- zhr_exons;

string <- as.character(obj[,1]);

exonLengths <- cbind(
do.call(rbind,strsplit(string,split="_"))[,1],
do.call(rbind,strsplit(do.call(rbind,strsplit(string,split="_"))[,2],split=","))[,1],
do.call(rbind,strsplit(do.call(rbind,strsplit(string,split="_"))[,2],split=","))[,2]
);

temp <- cbind(as.numeric(exonLengths[,2]),as.numeric(exonLengths[,3]))
length <- temp[,2]-temp[,1]+1;
depth <- obj[,2];
coverage <- (depth/length)*76;

hist(coverage,xlim=c(0,100),breaks=200);