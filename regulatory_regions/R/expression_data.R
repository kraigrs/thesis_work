library(vioplot);

zhrXz30 <- read.table("/Users/kraigrs/Wittkopp/regulatory_regions/data/zhrXz30_mRNA_meta_gene_output_121312_classified_C.txt",header=TRUE,sep="\t");
tsimbazazaXdroSec1 <- read.table("/Users/kraigrs/Wittkopp/regulatory_regions/data/tsimbazazaXdroSec1_mRNA_meta_gene_output_121312_classified_C.txt",header=TRUE,sep="\t");
zhrXtsimbazaza <- read.table("/Users/kraigrs/Wittkopp/regulatory_regions/data/zhrXtsimbazaza_mRNA_meta_gene_output_121312_classified_C.txt",header=TRUE,sep="\t");

# compare magnitude of expression changes across species comparisons for cis and trans

cis1 <- subset(zhrXz30, allcis == 1);
trans1 <- subset(zhrXz30, alltrans == 1);
cis2 <- subset(tsimbazazaXdroSec1, allcis == 1);
trans2 <- subset(tsimbazazaXdroSec1, alltrans == 1);
cis3 <- subset(zhrXtsimbazaza, allcis == 1);
trans3 <- subset(zhrXtsimbazaza, alltrans == 1);

#cis1 <- subset(zhrXz30, cis == 1);
#trans1 <- subset(zhrXz30, trans == 1);
#cis2 <- subset(tsimbazazaXdroSec1, cis == 1);
#trans2 <- subset(tsimbazazaXdroSec1, trans == 1);
#cis3 <- subset(zhrXtsimbazaza, cis == 1);
#trans3 <- subset(zhrXtsimbazaza, trans == 1);

v1 <- abs(log2(cis1$H1.s1/cis1$H1.s2));
v2 <- abs(log2(trans1$H1.s1/trans1$H1.s2));
v3 <- abs(log2(cis2$H1.s1/cis2$H1.s2));
v4 <- abs(log2(trans2$H1.s1/trans2$H1.s2));
v5 <- abs(log2(cis3$H1.s1/cis3$H1.s2));
v6 <- abs(log2(trans3$H1.s1/trans3$H1.s2));

v1 <- v1[!is.na(v1)];
v2 <- v2[!is.na(v2)];
v3 <- v3[!is.na(v3)];
v4 <- v4[!is.na(v4)];
v5 <- v5[!is.na(v5)];
v6 <- v6[!is.na(v6)];

vioplot(v1,v2,v3,v4,v5,v6,
names = rep(c("cis","trans"),3),
col="cyan",
ylim=c(0,10));

# compare magnitude of expression changes in categories within species comparison

data <- zhrXz30;

cis <- subset(data, allcis == 1);
trans <- subset(data, alltrans == 1);
compensatory <- subset(data, compensatory == 1);
cisplustrans <- subset(data, cisplustrans == 1);
cisXtrans <- subset(data, cisXtrans == 1);
conserved <- subset(data, conserved1 == 1);
ambiguous <- subset(data, ambiguous == 1);

v1 <- abs(log2(cis$H1.s1/cis$H1.s2));
v2 <- abs(log2(trans$H1.s1/trans$H1.s2));
v3 <- abs(log2(compensatory$H1.s1/compensatory$H1.s2));
v4 <- abs(log2(cisplustrans$H1.s1/cisplustrans$H1.s2));
v5 <- abs(log2(cisXtrans$H1.s1/cisXtrans$H1.s2));
v6 <- abs(log2(conserved$H1.s1/conserved$H1.s2));
v7 <- abs(log2(ambiguous$H1.s1/ambiguous$H1.s2));

v1 <- v1[!is.na(v1)];
v2 <- v2[!is.na(v2)];
v3 <- v3[!is.na(v3)];
v4 <- v4[!is.na(v4)];
v5 <- v5[!is.na(v5)];
v6 <- v6[!is.na(v6)];
v7 <- v7[!is.na(v7)];

vioplot(v1,v2,v3,v4,v5,v6,v7,
names = c("cis","trans","compensatory","cis+trans","cisXtrans","conserved","ambiguous"),
col="cyan",
ylim=c(0,10));

# compare magnitude expression to regulatory sequence divergence

data <- zhrXz30;
reg <- read.table("/Users/kraigrs/Wittkopp/regulatory_regions/references/zhr_z30/dm3_regulatory_regions_merged.zhr_z30_fsa.genes.txt",header=TRUE,sep="\t");
merged1 <- merge(reg,data,by.x="gene",by.y="gene");

data <- tsimbazazaXdroSec1;
reg <- read.table("/Users/kraigrs/Wittkopp/regulatory_regions/references/tsim_droSec1/dm3_regulatory_regions_merged.tsim_droSec1_fsa.genes.txt",header=TRUE,sep="\t");
merged2 <- merge(reg,data,by.x="gene",by.y="gene");

data <- zhrXtsimbazaza;
reg <- read.table("/Users/kraigrs/Wittkopp/regulatory_regions/references/zhr_tsim/dm3_regulatory_regions_merged.zhr_tsim_fsa.genes.txt",header=TRUE,sep="\t");
merged3 <- merge(reg,data,by.x="gene",by.y="gene");

par(mfrow=c(2,3))

plot(merged1$intronic_sites/(merged1$intronic_length/1000), abs(log2(merged1$H1.s1/merged1$H1.s2)),pch=19,cex=0.3,col=rgb(0,0,0,0.1),xlab="",ylab="",main="",xlim=c(0,20),ylim=c(0,2));
plot(merged2$intronic_sites/(merged2$intronic_length/1000), abs(log2(merged2$H1.s1/merged2$H1.s2)),pch=19,cex=0.3,col=rgb(0,0,0,0.1),xlab="",ylab="",main="",xlim=c(0,20),ylim=c(0,2));
plot(merged3$intronic_sites/(merged3$intronic_length/1000), abs(log2(merged3$H1.s1/merged3$H1.s2)),pch=19,cex=0.3,col=rgb(0,0,0,0.1),xlab="",ylab="",main="",xlim=c(0,20),ylim=c(0,2));

plot(merged1$intergenic_sites/(merged1$intergenic_length/1000), abs(log2(merged1$H1.s1/merged1$H1.s2)),pch=19,cex=0.3,col=rgb(0,0,0,0.1),xlab="",ylab="",main="",xlim=c(0,20),ylim=c(0,2));
plot(merged2$intergenic_sites/(merged2$intergenic_length/1000), abs(log2(merged2$H1.s1/merged2$H1.s2)),pch=19,cex=0.3,col=rgb(0,0,0,0.1),xlab="",ylab="",main="",xlim=c(0,20),ylim=c(0,2));
plot(merged3$intergenic_sites/(merged3$intergenic_length/1000), abs(log2(merged3$H1.s1/merged3$H1.s2)),pch=19,cex=0.3,col=rgb(0,0,0,0.1),xlab="",ylab="",main="",xlim=c(0,20),ylim=c(0,2));



lm( abs(log2(H1.s1/H1.s2)) ~ intronic_sites/(intronic_length/1000), data = merged1);
x <- merged1$intronic_sites/(merged1$intronic_length/1000);
y <- abs(log2(merged1$H1.s1/merged1$H1.s2));
linear <- lm(y~x);
summary(linear);

