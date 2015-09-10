# get expression values from mel_mel_data for simulations

zhr_exon <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/mRNA-Seq/zhr/zhr_v3.mosaik.exons.txt",sep="\t",header=TRUE);
z30_exon <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/mRNA-Seq/z30/z30_v3.mosaik.exons.txt",sep="\t",header=TRUE);

exon <- merge(zhr_exon,z30_exon,by.x="gene_exon",by.y="gene_exon");

SD_correct <- 16464075/21806797; # zhr / z30 --> z30 was sequenced at a greater depth

exprn <- cbind(exon,exon$zhr.x,round(exon$z30.y*SD_correct),round(exon$Both.x+exon$Both.y*SD_correct),round(exon$zhr.x+exon$z30.y*SD_correct+exon$Both.x+exon$Both.y*SD_correct));

total <- exprn[,c(1,15)];

write.table(total,file="/Users/kraigrs/Wittkopp/Simulations/zhr_z30_exons_expression.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE);
