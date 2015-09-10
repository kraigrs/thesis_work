############################################################################################
# script to look at sequence divergence
############################################################################################

annotation <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/genomes/dpse_r3.1.genes.bed",header=FALSE,sep="\t");
colnames(annotation) <- c("chromosome","start","stop","gene","empty","strand");

seq_div <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/genomes/dpse_r3.1.Dpse_TL.Dbog_Toro1.BQSR_genes.seq_div.txt",header=TRUE,sep="\t");

info <- merge(annotation,seq_div,by.x="gene",by.y="gene");

reliable <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/reliable_genes.txt",header=TRUE,sep="\t");

data_annotation <- merge(info,reliable,by.x="gene",by.y="gene");

# all genes

data_XL <- data_annotation[grep("^XL",data_annotation$chromosome,perl=TRUE),];
data_4 <- data_annotation[grep("^4",data_annotation$chromosome,perl=TRUE),];
data_3 <- data_annotation[grep("^3",data_annotation$chromosome,perl=TRUE),];
data_XR <- data_annotation[grep("^XR",data_annotation$chromosome,perl=TRUE),];
data_2 <- data_annotation[grep("^2",data_annotation$chromosome,perl=TRUE),];

boxplot((data_XL$SNPs/data_XL$coding)*100,at=0.25,xlim=c(0,2.5),ylim=c(0,2),xlab="",ylab="% sequence divergence",main="",outline=FALSE);
boxplot((data_4$SNPs/data_4$coding)*100,add=TRUE,at=0.75,xlim=c(0,2.5),ylim=c(0,2),xlab="",ylab="",main="",yaxt="n",outline=FALSE);
boxplot((data_3$SNPs/data_3$coding)*100,add=TRUE,at=1.25,xlim=c(0,2.5),ylim=c(0,2),xlab="",ylab="",main="",yaxt="n",outline=FALSE);
boxplot((data_XR$SNPs/data_XR$coding)*100,add=TRUE,at=1.75,xlim=c(0,2.5),ylim=c(0,2),xlab="",ylab="",main="",yaxt="n",outline=FALSE);
boxplot((data_2$SNPs/data_2$coding)*100,add=TRUE,at=2.25,xlim=c(0,2.5),ylim=c(0,2),xlab="",ylab="",main="",yaxt="n",outline=FALSE);
axis(side=1,at=c(0.25,0.75,1.25,1.75,2.25),labels=c("X","4","3","neo-X","2"));

# male-biased genes

#MBG <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/fold_change_carcass_MBG.txt",header=TRUE,sep="\t");
MBG <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/fold_change_gonads_MBG.txt",header=TRUE,sep="\t");
#MBG <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/testis_specific.txt",header=TRUE,sep="\t");
male <- merge(data_annotation,MBG,by.x="gene",by.y="gene");

male_XL <- male[grep("^XL",male$chromosome,perl=TRUE),];
male_4 <- male[grep("^4",male$chromosome,perl=TRUE),];
male_3 <- male[grep("^3",male$chromosome,perl=TRUE),];
male_XR <- male[grep("^XR",male$chromosome,perl=TRUE),];
male_2 <- male[grep("^2",male$chromosome,perl=TRUE),];

# female-biased genes

#FBG <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/fold_change_carcass_FBG.txt",header=TRUE,sep="\t");
FBG <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/fold_change_gonads_FBG.txt",header=TRUE,sep="\t");
#FBG <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/ovary_specific.txt",header=TRUE,sep="\t");
female <- merge(data_annotation,FBG,by.x="gene",by.y="gene");

female_XL <- female[grep("^XL",female$chromosome,perl=TRUE),];
female_4 <- female[grep("^4",female$chromosome,perl=TRUE),];
female_3 <- female[grep("^3",female$chromosome,perl=TRUE),];
female_XR <- female[grep("^XR",female$chromosome,perl=TRUE),];
female_2 <- female[grep("^2",female$chromosome,perl=TRUE),];

# non sex-biased genes

genes <- merge(male,female,all=TRUE)[,1];

index <- data_annotation$gene %in% genes;
data <- data_annotation[!index,];

data_XL <- data[grep("^XL",data$chromosome,perl=TRUE),];
data_4 <- data[grep("^4",data$chromosome,perl=TRUE),];
data_3 <- data[grep("^3",data$chromosome,perl=TRUE),];
data_XR <- data[grep("^XR",data$chromosome,perl=TRUE),];
data_2 <- data[grep("^2",data$chromosome,perl=TRUE),];

# make boxplots of sequence divergence

#boxplot((male_XL$SNPs/male_XL$coding)*100,at=0.1,xlim=c(0,9.2),ylim=c(0,3),xlab="",ylab="% sequence divergence",main="Carcass sex-biased genes",outline=FALSE,col="cornflowerblue");
boxplot((male_XL$SNPs/male_XL$coding)*100,at=0.1,xlim=c(0,9.2),ylim=c(0,3),xlab="",ylab="% sequence divergence",main="Gonad sex-biased genes",outline=FALSE,col="cornflowerblue");
#boxplot((male_XL$SNPs/male_XL$coding)*100,at=0.1,xlim=c(0,9.2),ylim=c(0,3),xlab="",ylab="% sequence divergence",main="Sex tissue-specific genes",outline=FALSE,col="cornflowerblue");
boxplot((data_XL$SNPs/data_XL$coding)*100,add=TRUE,at=0.6,xlim=c(0,9.2),ylim=c(0,3),xlab="",ylab="",main="",yaxt="n",outline=FALSE);
boxplot((female_XL$SNPs/female_XL$coding)*100,add=TRUE,at=1.1,xlim=c(0,9.2),ylim=c(0,3),xlab="",ylab="",main="",yaxt="n",outline=FALSE,col="plum1");

boxplot((male_4$SNPs/male_4$coding)*100,add=TRUE,at=2.1,xlim=c(0,9.2),ylim=c(0,3),xlab="",ylab="",main="",yaxt="n",outline=FALSE,col="cornflowerblue");
boxplot((data_4$SNPs/data_4$coding)*100,add=TRUE,at=2.6,xlim=c(0,9.2),ylim=c(0,3),xlab="",ylab="",main="",yaxt="n",outline=FALSE);
boxplot((female_4$SNPs/female_4$coding)*100,add=TRUE,at=3.1,xlim=c(0,9.2),ylim=c(0,3),xlab="",ylab="",main="",yaxt="n",outline=FALSE,col="plum1");

boxplot((male_3$SNPs/male_3$coding)*100,add=TRUE,at=4.1,xlim=c(0,9.2),ylim=c(0,3),xlab="",ylab="",main="",yaxt="n",outline=FALSE,col="cornflowerblue");
boxplot((data_3$SNPs/data_3$coding)*100,add=TRUE,at=4.6,xlim=c(0,9.2),ylim=c(0,3),xlab="",ylab="",main="",yaxt="n",outline=FALSE);
boxplot((female_3$SNPs/female_3$coding)*100,add=TRUE,at=5.1,xlim=c(0,9.2),ylim=c(0,3),xlab="",ylab="",main="",yaxt="n",outline=FALSE,col="plum1");

boxplot((male_XR$SNPs/male_XR$coding)*100,add=TRUE,at=6.1,xlim=c(0,9.2),ylim=c(0,3),xlab="",ylab="",main="",yaxt="n",outline=FALSE,col="cornflowerblue");
boxplot((data_XR$SNPs/data_XR$coding)*100,add=TRUE,at=6.6,xlim=c(0,9.2),ylim=c(0,3),xlab="",ylab="",main="",yaxt="n",outline=FALSE);
boxplot((female_XR$SNPs/female_XR$coding)*100,add=TRUE,at=7.1,xlim=c(0,9.2),ylim=c(0,3),xlab="",ylab="",main="",yaxt="n",outline=FALSE,col="plum1");

boxplot((male_2$SNPs/male_2$coding)*100,add=TRUE,at=8.1,xlim=c(0,9.2),ylim=c(0,3),xlab="",ylab="",main="",yaxt="n",outline=FALSE,col="cornflowerblue");
boxplot((data_2$SNPs/data_2$coding)*100,add=TRUE,at=8.6,xlim=c(0,9.2),ylim=c(0,3),xlab="",ylab="",main="",yaxt="n",outline=FALSE);
boxplot((female_2$SNPs/female_2$coding)*100,add=TRUE,at=9.1,xlim=c(0,9.2),ylim=c(0,3),xlab="",ylab="",main="",yaxt="n",outline=FALSE,col="plum1");

axis(side=1,at=c(0.6,2.6,4.6,6.6,8.6),labels=c("X","4","3","neo-X","2"));

legend("topright",legend=c("male","female"),fill=c("cornflowerblue","plum1"));

