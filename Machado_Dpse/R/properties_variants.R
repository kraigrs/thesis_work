Dpse_TL_vars <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/genomes/Dpse_TL.GATK_variants.genes.txt",header=TRUE,sep="\t");
Dpse_TL_expr <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/Dpse_TL.genes.txt",header=TRUE,sep="\t");

Dbog_Toro1_vars <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/genomes/Dbog_Toro1.GATK_variants.genes.txt",header=TRUE,sep="\t");
Dbog_Toro1_expr <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/Dbog_Toro1.genes.txt",header=TRUE,sep="\t");

Dpse_TL <- merge(Dpse_TL_vars,Dpse_TL_expr,by.x="gene",by.y="gene");
Dbog_Toro1 <- merge(Dbog_Toro1_vars,Dbog_Toro1_expr,by.x="gene",by.y="gene");

data <- merge(Dpse_TL,Dbog_Toro1,by.x="gene",by.y="gene");

annotation <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/genomes/dpse_r3.1.genes.bed",header=FALSE,sep="\t");
colnames(annotation) <- c("chromosome","start","stop","gene","empty","strand");

seq_div <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/genomes/dpse_r3.1.Dpse_TL.Dbog_Toro1.BQSR_genes.seq_div.txt",header=TRUE,sep="\t");

info <- merge(annotation,seq_div,by.x="gene",by.y="gene");

data_annotation <- merge(data,info,by.x="gene",by.y="gene");

temp <- subset(data_annotation,Dpse_TL.x+Dbog_Toro1.y>=20);
data_annotation <- temp;

plot(data_annotation$Dpse_TL.x/(data_annotation$Dpse_TL.x+data_annotation$Dbog_Toro1.x),
     data_annotation$Dbog_Toro1.y/(data_annotation$Dpse_TL.y+data_annotation$Dbog_Toro1.y),
     xlab="Dpse gDNA correct mapping",ylab="Dbog gDNA correct mapping",
     pch=16,cex=0.75,col=rgb(0,0,0,0.5));

data_XL <- data_annotation[grep("^XL",data_annotation$chromosome,perl=TRUE),];
data_4 <- data_annotation[grep("^4",data_annotation$chromosome,perl=TRUE),];
data_3 <- data_annotation[grep("^3",data_annotation$chromosome,perl=TRUE),];
data_XR <- data_annotation[grep("^XR",data_annotation$chromosome,perl=TRUE),];
data_2 <- data_annotation[grep("^2",data_annotation$chromosome,perl=TRUE),];



######## produce heatmap-style contour plot ########

library(ggplot2);
library(RColorBrewer);

k <- 11;
my.cols <- rev(brewer.pal(k, "RdYlBu"));
n <- nrow(subset(data_annotation,Dpse_TL.x/(Dpse_TL.x+Dbog_Toro1.x)<0.8 | Dbog_Toro1.y/(Dpse_TL.y+Dbog_Toro1.y)<0.8));
n <- nrow(data_annotation);

pdf(file="/Users/kraigrs/Desktop/gDNA_correct_mapping.pdf",height=5,width=10);

par(mfrow=c(1,2));

plot(data_annotation$Dpse_TL.x/(data_annotation$Dpse_TL.x+data_annotation$Dbog_Toro1.x),
     data_annotation$Dbog_Toro1.y/(data_annotation$Dpse_TL.y+data_annotation$Dbog_Toro1.y),
     xlab="Dpse gDNA correct mapping",ylab="Dbog gDNA correct mapping",
     pch=16,cex=0.75,col=rgb(0,0,0,0.5));

smoothScatter(data_annotation$Dpse_TL.x/(data_annotation$Dpse_TL.x+data_annotation$Dbog_Toro1.x),data_annotation$Dbog_Toro1.y/(data_annotation$Dpse_TL.y+data_annotation$Dbog_Toro1.y),pch=19,cex=0.8,col=rgb(0,0,0,0.5),nrpoints=0, colramp=colorRampPalette(my.cols),xlab="Dpse gDNA correct mapping",ylab="Dbog gDNA correct mapping");

dev.off();

####################################################

data_mismap <- subset(data_annotation,Dpse_TL.x/(Dpse_TL.x+Dbog_Toro1.x)<0.9 | Dbog_Toro1.y/(Dpse_TL.y+Dbog_Toro1.y)<0.9);
data_perfect_in_one <- subset(data_annotation,Dpse_TL.x/(Dpse_TL.x+Dbog_Toro1.x)==1 | Dbog_Toro1.y/(Dpse_TL.y+Dbog_Toro1.y)==1);
data_perfect <- subset(data_annotation,Dpse_TL.x/(Dpse_TL.x+Dbog_Toro1.x)==1 & Dbog_Toro1.y/(Dpse_TL.y+Dbog_Toro1.y)==1);
data_filter <- subset(data_annotation,Dpse_TL.x/(Dpse_TL.x+Dbog_Toro1.x)>=0.9 & Dbog_Toro1.y/(Dpse_TL.y+Dbog_Toro1.y)>=0.9);

plot(data_mismap$Dpse_TL.x/(data_mismap$Dpse_TL.x+data_mismap$Dbog_Toro1.x),
     1-data_mismap$Dbog_Toro1.y/(data_mismap$Dpse_TL.y+data_mismap$Dbog_Toro1.y),
     xlab="Dpse gDNA mismapping",ylab="Dbog gDNA mismapping",
     pch=16,cex=0.75,col=rgb(0,0,0,0.5));

data_XL <- data_annotation[grep("^XL",data_annotation$chromosome,perl=TRUE),];
data_4 <- data_annotation[grep("^4",data_annotation$chromosome,perl=TRUE),];
data_3 <- data_annotation[grep("^3",data_annotation$chromosome,perl=TRUE),];
data_XR <- data_annotation[grep("^XR",data_annotation$chromosome,perl=TRUE),];
data_2 <- data_annotation[grep("^2",data_annotation$chromosome,perl=TRUE),];

boxplot((data_XL$SNPs/data_XL$coding)*100,at=0.5,xlim=c(0,3),ylim=c(0,3),xlab="",ylab="",main="",yaxt="n",pch=16);
boxplot((data_4$SNPs/data_4$coding)*100,add=TRUE,at=1,xlim=c(0,3),ylim=c(0,3),xlab="",ylab="",main="",yaxt="n",pch=16);
boxplot((data_3$SNPs/data_3$coding)*100,add=TRUE,at=1.5,xlim=c(0,3),ylim=c(0,3),xlab="",ylab="",main="",yaxt="n",pch=16);
boxplot((data_XR$SNPs/data_XR$coding)*100,add=TRUE,at=2,xlim=c(0,3),ylim=c(0,3),xlab="",ylab="",main="",yaxt="n",pch=16);
boxplot((data_2$SNPs/data_2$coding)*100,add=TRUE,at=2.5,xlim=c(0,3),ylim=c(0,3),xlab="",ylab="",main="",yaxt="n",pch=16);
axis(side=1,at=c(0.5,1,1.5,2,2.5),labels=c("X","4","3","neo-X","2"));

mean((data_XL$SNPs/data_XL$coding)*100);
mean((data_4$SNPs/data_4$coding)*100);
mean((data_3$SNPs/data_3$coding)*100);
mean((data_XR$SNPs/data_XR$coding)*100);
mean((data_2$SNPs/data_2$coding)*100);


##################
# homozygous variants

par(mfrow=c(1,2));

boxplot(data_XL$pass_hom.x/data_XL$coding,at=0.5,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",pch=16,outline=FALSE);
axis(side=1,at=c(0.75,2.25,3.75,5.25,6.75),labels=c("X","4","3","neo-X","2"));
legend("topleft",legend=c("Dpse","Dbog"),fill=c("white","gray"));
boxplot(data_XL$pass_hom.y/data_XL$coding,add=TRUE,at=1,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",col="gray",pch=16,outline=FALSE);

boxplot(data_4$pass_hom.x/data_4$coding,add=TRUE,at=2,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",pch=16,outline=FALSE);
boxplot(data_4$pass_hom.y/data_4$coding,add=TRUE,at=2.5,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",col="gray",pch=16,outline=FALSE);

boxplot(data_3$pass_hom.x/data_3$coding,add=TRUE,at=3.5,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",pch=16,outline=FALSE);
boxplot(data_3$pass_hom.y/data_3$coding,add=TRUE,at=4,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",col="gray",pch=16,outline=FALSE);

boxplot(data_XR$pass_hom.x/data_XR$coding,add=TRUE,at=5,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",pch=16,outline=FALSE);
boxplot(data_XR$pass_hom.y/data_XR$coding,add=TRUE,at=5.5,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",col="gray",pch=16,outline=FALSE);

boxplot(data_2$pass_hom.x/data_2$coding,add=TRUE,at=6.5,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",pch=16,outline=FALSE);
boxplot(data_2$pass_hom.y/data_2$coding,add=TRUE,at=7,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="# variants/coding bases",main="Homozygous variants passing filters",col="gray",pch=16,outline=FALSE);


boxplot(data_XL$fail_hom.x/data_XL$coding,at=0.5,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",pch=16,outline=FALSE);
axis(side=1,at=c(0.75,2.25,3.75,5.25,6.75),labels=c("X","4","3","neo-X","2"));
legend("topleft",legend=c("Dpse","Dbog"),fill=c("white","gray"));
boxplot(data_XL$fail_hom.y/data_XL$coding,add=TRUE,at=1,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",col="gray",pch=16,outline=FALSE);

boxplot(data_4$fail_hom.x/data_4$coding,add=TRUE,at=2,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",pch=16,outline=FALSE);
boxplot(data_4$fail_hom.y/data_4$coding,add=TRUE,at=2.5,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",col="gray",pch=16,outline=FALSE);

boxplot(data_3$fail_hom.x/data_3$coding,add=TRUE,at=3.5,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",pch=16,outline=FALSE);
boxplot(data_3$fail_hom.y/data_3$coding,add=TRUE,at=4,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",col="gray",pch=16,outline=FALSE);

boxplot(data_XR$fail_hom.x/data_XR$coding,add=TRUE,at=5,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",pch=16,outline=FALSE);
boxplot(data_XR$fail_hom.y/data_XR$coding,add=TRUE,at=5.5,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",col="gray",pch=16,outline=FALSE);

boxplot(data_2$fail_hom.x/data_2$coding,add=TRUE,at=6.5,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",pch=16,outline=FALSE);
boxplot(data_2$fail_hom.y/data_2$coding,add=TRUE,at=7,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="Homozygous variants failing filters",col="gray",pch=16,outline=FALSE);


######################
# heterozygous variants

boxplot(data_XL$pass_het.x/data_XL$coding,at=0.5,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",pch=16,outline=FALSE);
axis(side=1,at=c(0.75,2.25,3.75,5.25,6.75),labels=c("X","4","3","neo-X","2"));
boxplot(data_XL$pass_het.y/data_XL$coding,add=TRUE,at=1,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",col="gray",pch=16,outline=FALSE);

boxplot(data_4$pass_het.x/data_4$coding,add=TRUE,at=2,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",pch=16,outline=FALSE);
boxplot(data_4$pass_het.y/data_4$coding,add=TRUE,at=2.5,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",col="gray",pch=16,outline=FALSE);

boxplot(data_3$pass_het.x/data_3$coding,add=TRUE,at=3.5,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",pch=16,outline=FALSE);
boxplot(data_3$pass_het.y/data_3$coding,add=TRUE,at=4,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",col="gray",pch=16,outline=FALSE);

boxplot(data_XR$pass_het.x/data_XR$coding,add=TRUE,at=5,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",pch=16,outline=FALSE);
boxplot(data_XR$pass_het.y/data_XR$coding,add=TRUE,at=5.5,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",col="gray",pch=16,outline=FALSE);

boxplot(data_2$pass_het.x/data_2$coding,add=TRUE,at=6.5,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",pch=16,outline=FALSE);
boxplot(data_2$pass_het.y/data_2$coding,add=TRUE,at=7,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="# sites/coding bases",main="Heterozygous variants passing filters",col="gray",pch=16,outline=FALSE);


boxplot(data_XL$fail_het.x/data_XL$coding,at=0.5,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",pch=16,outline=FALSE);
axis(side=1,at=c(0.75,2.25,3.75,5.25,6.75),labels=c("X","4","3","neo-X","2"));
boxplot(data_XL$fail_het.y/data_XL$coding,add=TRUE,at=1,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",col="gray",pch=16,outline=FALSE);

boxplot(data_4$fail_het.x/data_4$coding,add=TRUE,at=2,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",pch=16,outline=FALSE);
boxplot(data_4$fail_het.y/data_4$coding,add=TRUE,at=2.5,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",col="gray",pch=16,outline=FALSE);

boxplot(data_3$fail_het.x/data_3$coding,add=TRUE,at=3.5,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",pch=16,outline=FALSE);
boxplot(data_3$fail_het.y/data_3$coding,add=TRUE,at=4,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",col="gray",pch=16,outline=FALSE);

boxplot(data_XR$fail_het.x/data_XR$coding,add=TRUE,at=5,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",pch=16,outline=FALSE);
boxplot(data_XR$fail_het.y/data_XR$coding,add=TRUE,at=5.5,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",yaxt="n",col="gray",pch=16,outline=FALSE);

boxplot(data_2$fail_het.x/data_2$coding,add=TRUE,at=6.5,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="",pch=16,outline=FALSE);
boxplot(data_2$fail_het.y/data_2$coding,add=TRUE,at=7,ylim=c(0,0.02),xlim=c(0,7.5),xlab="",ylab="",main="Heterozygous variants failing filters",col="gray",pch=16,outline=FALSE);


#################################################
# combine mismapping and variants failing filters

map_XL <- subset(data_XL,Dpse_TL.x/(Dpse_TL.x+Dbog_Toro1.x) >= 0.9 & Dbog_Toro1.y/(Dpse_TL.y+Dbog_Toro1.y) >= 0.9);
mismap_XL <- subset(data_XL,Dpse_TL.x/(Dpse_TL.x+Dbog_Toro1.x) < 0.9 | Dbog_Toro1.y/(Dpse_TL.y+Dbog_Toro1.y) < 0.9);



nrow(subset(data_XL,Dpse_TL.x/(Dpse_TL.x+Dbog_Toro1.x) >= 0.9 & Dbog_Toro1.y/(Dpse_TL.y+Dbog_Toro1.y) >= 0.9))/nrow(data_XL);
nrow(subset(data_4,Dpse_TL.x/(Dpse_TL.x+Dbog_Toro1.x) >= 0.9 & Dbog_Toro1.y/(Dpse_TL.y+Dbog_Toro1.y) >= 0.9))/nrow(data_4);
nrow(subset(data_3,Dpse_TL.x/(Dpse_TL.x+Dbog_Toro1.x) >= 0.9 & Dbog_Toro1.y/(Dpse_TL.y+Dbog_Toro1.y) >= 0.9))/nrow(data_3);
nrow(subset(data_XR,Dpse_TL.x/(Dpse_TL.x+Dbog_Toro1.x) >= 0.9 & Dbog_Toro1.y/(Dpse_TL.y+Dbog_Toro1.y) >= 0.9))/nrow(data_XR);
nrow(subset(data_2,Dpse_TL.x/(Dpse_TL.x+Dbog_Toro1.x) >= 0.9 & Dbog_Toro1.y/(Dpse_TL.y+Dbog_Toro1.y) >= 0.9))/nrow(data_2);


cor(data_XL$fail_hom.x/data_XL$coding,1-data_XL$Dpse_TL.x/(data_XL$Dpse_TL.x+data_XL$Dbog_Toro1.x),use="pairwise.complete.obs",method="spearman");
cor(data_4$fail_hom.x/data_4$coding,1-data_4$Dpse_TL.x/(data_4$Dpse_TL.x+data_4$Dbog_Toro1.x),use="pairwise.complete.obs",method="spearman");
cor(data_3$fail_hom.x/data_3$coding,1-data_3$Dpse_TL.x/(data_3$Dpse_TL.x+data_3$Dbog_Toro1.x),use="pairwise.complete.obs",method="spearman");
cor(data_XR$fail_hom.x/data_XR$coding,1-data_XR$Dpse_TL.x/(data_XR$Dpse_TL.x+data_XR$Dbog_Toro1.x),use="pairwise.complete.obs",method="spearman");
cor(data_2$fail_hom.x/data_2$coding,1-data_2$Dpse_TL.x/(data_2$Dpse_TL.x+data_2$Dbog_Toro1.x),use="pairwise.complete.obs",method="spearman");

plot(data_annotation$fail_hom.x/data_annotation$coding,1-data_annotation$Dpse_TL.x/(data_annotation$Dpse_TL.x+data_annotation$Dbog_Toro1.x),xlab="# failed variants/coding bases",ylab="proportion mismapping",pch=16,col=rgb(0,0,0,0.5),cex=0.5);








