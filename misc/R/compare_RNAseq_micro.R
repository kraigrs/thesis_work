genes <- read.table("/Users/kraigrs/Wittkopp/McManus/constitutive_genes.txt",sep="\t",na.strings=c("-"," "),fill=TRUE);

genes <- subset(genes,!is.na(genes[,3]));
genes <- subset(genes,!is.na(genes[,4]));
lengths <- genes[,3]-genes[,4];
genes <- cbind(genes,lengths);
genes <- genes[,c(1,ncol(genes))];

zhr_rna_seq <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/mRNA-Seq/zhr/zhr.mosaik.genes.txt",sep="\t",header=TRUE);
z30_rna_seq <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/mRNA-Seq/z30/z30.mosaik.genes.txt",sep="\t",header=TRUE);

total_zhr <- zhr_rna_seq$zhr + zhr_rna_seq$z30 + zhr_rna_seq$Both;
zhr_rna_seq <- cbind(zhr_rna_seq,total_zhr);
zhr_rna_seq <- zhr_rna_seq[,c(1,ncol(zhr_rna_seq))];

total_z30 <- z30_rna_seq$zhr + z30_rna_seq$z30 + z30_rna_seq$Both;
z30_rna_seq <- cbind(z30_rna_seq,total_z30);
z30_rna_seq <- z30_rna_seq[,c(1,ncol(z30_rna_seq))];

rna_seq <- merge(zhr_rna_seq,z30_rna_seq,by.x="gene",by.y="gene"); 

micro <- read.table("/Users/kraigrs/Wittkopp/DGRP/The40TranscriptomeExpressionData.txt",sep="\t",skip=2);
temp <- cbind(micro,2^(micro[,2:ncol(micro)]));
average <- rowMeans(temp[,82:ncol(micro)]);
micro <- cbind(temp,average);
#micro <- micro[,1:41];

#average <- rowMeans(micro[,2:ncol(micro)]);
micro <- cbind(micro,average);
micro <- micro[,c(1,ncol(micro))];

# merge all data

map <- read.table("/Users/kraigrs/Wittkopp/DGRP/DGRP_probe2gene_map.txt",header=TRUE,sep="\t");

m1 <- merge(rna_seq,genes,by.x="gene",by.y="V1");
m2 <- merge(m1,map,by.x="gene",by.y="gene");

data <- merge(m2,micro,by.x="probe",by.y="V1"); # 8,978 rows
data <- subset(data,data$total_zhr > 0); 
data <- subset(data,data$total_z30 > 0); # 7,876 rows

plot(log2(data$total_z30),log2(data$total_zhr),xlab="log2(total z30)",ylab="log2(total zhr)",main="",xlim=c(0,16.5),ylim=c(0,16.5));

# RPKM
RPKM_zhr <- (data$total_zhr*1e9)/(data$lengths*13156445);
RPKM_z30 <- (data$total_z30*1e9)/(data$lengths*18536664);

par(mfrow=c(1,2),oma = c( 0, 0, 3, 0 ) );

#DGRP microarray results from 40 female lines of Dmel\ncompared to RNA-seq results from Dmel strains

plot(log2(data$total_zhr),data$average,xlab="Illumina (log2) counts",ylab="Affymetrix intensities",main="zhr",xlim=c(0,16.5),ylim=c(8.9,16.5),pch=19,col=rgb(0,0,0,0.3),cex=0.4);
plot(log2(data$total_z30),data$average,xlab="Illumina (log2) counts",ylab="Affymetrix intensities",main="z30",xlim=c(0,16.5),ylim=c(8.9,16.5),pch=19,col=rgb(0,0,0,0.3),cex=0.4);

title("DGRP microarray results from 40 female lines of Dmel\ncompared to RNA-seq results from Dmel strains", outer = TRUE );

# separate into tertiles

firstTert <- subset(data,data$average <= quantile(data$average,probs=1/3));
secondTert <- subset(data,data$average > quantile(data$average,probs=1/3) & data$average <= quantile(data$average,probs=2/3));
thirdTert <- subset(data,data$average > quantile(data$average,probs=2/3));

level <- rep(1,nrow(firstTert));
firstTert <- cbind(firstTert,level);

level <- rep(5,nrow(secondTert));
secondTert <- cbind(secondTert,level);

level <- rep(10,nrow(thirdTert));
thirdTert <- cbind(thirdTert,level);

data <- rbind(firstTert,secondTert,thirdTert);

#compare number of transcripts mapped to number of total possible transcripts

#compare exon length with total expression
zhr_exon <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/mRNA-Seq/zhr/zhr_v3.mosaik.exons.txt",sep="\t",header=TRUE);
start <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/mRNA-Seq/zhr/start");
stop <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/mRNA-Seq/zhr/stop");
length <- stop-start;
zhr_exon <- cbind(zhr_exon,length);
total <- zhr_exon$zhr+zhr_exon$Both;
zhr_exon <- cbind(zhr_exon,total);

z30_exon <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/mRNA-Seq/z30/z30_v3.mosaik.exons.txt",sep="\t",header=TRUE);
start <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/mRNA-Seq/z30/start");
stop <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/mRNA-Seq/z30/stop");
length <- stop-start;
z30_exon <- cbind(z30_exon,length);
total <- z30_exon$z30+z30_exon$Both;
z30_exon <- cbind(z30_exon,total);

exon <- merge(zhr_exon,z30_exon,by.x="gene_exon",by.y="gene_exon");
# average of average adjusted zhr and z30 values, floored
avg <- floor(((exon$Adjzhr.x+exon$Adjzhr.y)/2 + (exon$Adjz30.x+exon$Adjz30.y)/2)/2);
exon <- cbind(exon,avg);
temp <- exon[,c(1,ncol(exon))];
write.table(temp,file="/Users/kraigrs/Wittkopp/mel_mel_data/mRNA-Seq/zhr_z30_exons_expression.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE);

par(mfrow=c(1,2));
plot(log2(zhr_exon$total),zhr_exon$V1,xlab="log2(transcript abundance)",ylab="exon length (bp)",main="zhr");
plot(log2(z30_exon$total),z30_exon$V1,xlab="log2(transcript abundance)",ylab="exon length (bp)",main="z30");

#get expression values for ASE simulation
zhrXz30 <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/mRNA-Seq/zhrXz30/zhrXz30.mosaik.exons.txt",sep="\t",header=TRUE);
z30Xzhr <- read.table("/Users/kraigrs/Wittkopp/mel_mel_data/mRNA-Seq/z30Xzhr/z30Xzhr.mosaik.exons.txt",sep="\t",header=TRUE);

data <- merge(zhrXz30,z30Xzhr,by.x="gene_exon",by.y="gene_exon");
zhr_avg <- floor((data$Adjzhr.x+data$Adjzhr.y)/2);
z30_avg <- floor((data$Adjz30.x+data$Adjz30.y)/2);
data <- cbind(data,zhr_avg,z30_avg);

m1 <- subset(data,data$zhr_avg + data$z30_avg >= 20);
m2 <- m1[,c(1,ncol(m1)-1,ncol(m1))];

write.table(m2,file="/Users/kraigrs/Wittkopp/mel_mel_data/mRNA-Seq/hybrid_exons_expression.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE);





