data1 <- read.table("/Users/kraigrs/Wittkopp/Graze/SRR389095.mate1.berlin-updated-exonic-regions.bowtie_v1_m1.exons.txt",header=TRUE,sep="\t");

data2 <- read.table("/Users/kraigrs/Wittkopp/Graze/SRR389095.mate1.berlin-updated-exonic-regions.bowtie_v2_m1.exons.txt",header=TRUE,sep="\t");

data3 <- read.table("/Users/kraigrs/Wittkopp/Graze/SRR389095.mate1.berlin-updated-exonic-regions.bowtie_v3_m1.exons.txt",header=TRUE,sep="\t");
data3_biased <- read.table("/Users/kraigrs/Wittkopp/Graze/SRR389095.mate1.berlin-updated-exonic-regions.bowtie_v3_m1.unfiltered.exons.txt",header=TRUE,sep="\t");

data0 <- read.table("/Users/kraigrs/Wittkopp/Graze/SRR389095.mate1.berlin_c1674_fsa_masked.bowtie_v0_m1.exons.txt",header=TRUE,sep="\t");

par(mfrow=c(2,2));
hist(data1$berlin/(data1$berlin+data1$c1674));
hist(data2$berlin/(data2$berlin+data2$c1674));
hist(data3$berlin/(data3$berlin+data3$c1674));
hist(data0$berlin/(data0$berlin+data0$c1674));

# AS >= 20

data1_complete <- subset(data1,berlin + c1674 >= 20);
data2_complete <- subset(data2,berlin + c1674 >= 20);
data3_complete <- subset(data3,berlin + c1674 >= 20);
data3_complete_biased <- subset(data3_biased,berlin + c1674 >= 20);
data0_complete <- subset(data0,berlin + c1674 >= 20);

d1 <- density(data1_complete$berlin/(data1_complete$berlin+data1_complete$c1674));
d2 <- density(data2_complete$berlin/(data2_complete$berlin+data2_complete$c1674));
d3 <- density(data3_complete$berlin/(data3_complete$berlin+data3_complete$c1674));
d0 <- density(data0_complete$berlin/(data0_complete$berlin+data0_complete$c1674));

par(mfrow=c(2,2));

plot(d1,xlim=c(0,1),ylim=c(0,8),main="single reference (1mm)");
plot(d2,xlim=c(0,1),ylim=c(0,8),main="single reference (2mm)");
plot(d3,xlim=c(0,1),ylim=c(0,8),main="single reference (3mm)");
plot(d0,xlim=c(0,1),ylim=c(0,8),main="multiple reference (0mm)");

pvals <- NULL;
for(i in 1:nrow(data1_complete))
{
	test <- binom.test(data1_complete$berlin[i],data1_complete$berlin[i]+data1_complete$c1674[i],alternative="two.sided");
	pvals <- rbind(pvals,test$p.value);
}
qvals <- p.adjust(pvals);
data1_tests <- cbind(data1_complete,pvals,qvals);
sum(data1_tests$qvals >= 0.05)/nrow(data1_tests);

pvals <- NULL;
for(i in 1:nrow(data2_complete))
{
	test <- binom.test(data2_complete$berlin[i],data2_complete$berlin[i]+data2_complete$c1674[i],alternative="two.sided");
	pvals <- rbind(pvals,test$p.value);
}
qvals <- p.adjust(pvals);
data2_tests <- cbind(data2_complete,pvals,qvals);
sum(data2_tests$qvals >= 0.05)/nrow(data2_tests);

pvals <- NULL;
for(i in 1:nrow(data3_complete))
{
	test <- binom.test(data3_complete$berlin[i],data3_complete$berlin[i]+data3_complete$c1674[i],alternative="two.sided");
	pvals <- rbind(pvals,test$p.value);
}
qvals <- p.adjust(pvals);
data3_tests <- cbind(data3_complete,pvals,qvals);
sum(data3_tests$qvals >= 0.05)/nrow(data3_tests);

pvals <- NULL;
for(i in 1:nrow(data3_complete_biased))
{
	test <- binom.test(data3_complete_biased$berlin[i],data3_complete_biased$berlin[i]+data3_complete_biased$c1674[i],alternative="two.sided");
	pvals <- rbind(pvals,test$p.value);
}
qvals <- p.adjust(pvals);
data3_tests_biased <- cbind(data3_complete_biased,pvals,qvals);
sum(data3_tests_biased$qvals >= 0.05)/nrow(data3_tests_biased);

pvals <- NULL;
for(i in 1:nrow(data0_complete))
{
	test <- binom.test(data0_complete$berlin[i],data0_complete$berlin[i]+data0_complete$c1674[i],alternative="two.sided");
	pvals <- rbind(pvals,test$p.value);
}
qvals <- p.adjust(pvals);
data0_tests <- cbind(data0_complete,pvals,qvals);
sum(data0_tests$qvals >= 0.05)/nrow(data0_tests);


compare <- merge(data0_tests,data3_tests_biased,by.x="exon",by.y="exon");
cor(compare$berlin.x/(compare$berlin.x+compare$c1674.x) , compare$berlin.y/(compare$berlin.y+compare$c1674.y))^2;
summary(compare$berlin.x+compare$c1674.x); 
summary(compare$berlin.y+compare$c1674.y);


nrow(compare);
s_ns <- subset(compare,qvals.x < 0.05 & qvals.y >= 0.05);
ns_s <- subset(compare,qvals.x >= 0.05 & qvals.y < 0.05);
tot <- nrow(s_ns)+nrow(ns_s);
a <- nrow(subset(s_ns,ref_allele/(ref_allele+alt_allele) > berlin_ref_allele/(berlin_ref_allele+c1674_alt_allele)));
b <- nrow(subset(ns_s,ref_allele/(ref_allele+alt_allele) > berlin_ref_allele/(berlin_ref_allele+c1674_alt_allele)));
(a+b)/tot;

#pdf(file="/Users/kraigrs/Desktop/new_plots/qqplot_real_data_qvals_3mm_unfiltered.pdf");
#qqplot(-log10(compare$qvals.x),-log10(compare$qvals.y),xlab="-log10(qvals) parental genomes",ylab="-log10(qvals) single genome",main="");
#abline(a=0,b=1,col="red",lty=2);
#dev.off();

pdf(file="/Users/kraigrs/Desktop/plots/compare_real_exons_3mm_biased.pdf");

compare <- merge(data0_tests,data3_tests_biased,by.x="exon",by.y="exon");
nrow(compare);

plot( compare$berlin.x/(compare$berlin.x+compare$c1674.x) , compare$berlin.y/(compare$berlin.y+compare$c1674.y) ,xlim=c(0,1),ylim=c(0,1),pch=19,cex=0.5,col=rgb(0,0,0,0.05),xlab="",ylab="",main="");

dev.off();



pdf(file="/Users/kraigrs/Desktop/plots/individual/pie_real_exons_3mm_biased.pdf");

compare <- merge(data0_tests,data3_tests_biased,by.x="exon",by.y="exon");

n_n <- subset(compare,qvals.x >= 0.05 & qvals.y >= 0.05);
s_ns <- subset(compare,qvals.x < 0.05 & qvals.y >= 0.05);
ns_s <- subset(compare,qvals.x >= 0.05 & qvals.y < 0.05);
s_s <- subset(compare,qvals.x < 0.05 & qvals.y < 0.05);

pie(c(nrow(n_n),nrow(s_ns),nrow(ns_s),nrow(s_s)),col=c("grey","red","blue","purple"),labels="");

dev.off();



ns_ns <- subset(compare,qvals.x >= 0.05 & qvals.y >= 0.05);
nrow(ns_ns);
points(ns_ns$berlin.x/(ns_ns$berlin.x+ns_ns$c1674.x),ns_ns$berlin.y/(ns_ns$berlin.y+ns_ns$c1674.y),pch=19,cex=0.5,col=rgb(0,0,0,0.05));

s_ns <- subset(compare,qvals.x < 0.05 & qvals.y >= 0.05);
nrow(s_ns);
points(s_ns$berlin.x/(s_ns$berlin.x+s_ns$c1674.x),s_ns$berlin.y/(s_ns$berlin.y +s_ns$c1674.y),pch=19,cex=0.5,col=rgb(1,0,0,0.6));

ns_s <- subset(compare,qvals.x >= 0.05 & qvals.y < 0.05);
nrow(ns_s);
points(ns_s$berlin.x/(ns_s$berlin.x+ns_s$c1674.x),ns_s$berlin.y/(ns_s$berlin.y +ns_s$c1674.y),pch=19,cex=0.5,col=rgb(0,0,1,0.6));

s_s <- subset(compare,qvals.x < 0.05 & qvals.y < 0.05);
nrow(s_s);
points(s_s$berlin.x/(s_s$berlin.x +s_s$c1674.x),s_s$berlin.y/(s_s$berlin.y +s_s$c1674.y),pch=19,cex=0.5,col=rgb(1,0,1,0.6));

dev.off();



# write to file so dont have to repeat tests

write.table(data1_tests,file="/Users/kraigrs/Wittkopp/Graze/SRR389095.mate1.berlin-updated-exonic-regions.bowtie_v1_m1.exons.binom_tests.txt",quote=FALSE,row.names=FALSE,sep="\t");

write.table(data2_tests,file="/Users/kraigrs/Wittkopp/Graze/SRR389095.mate1.berlin-updated-exonic-regions.bowtie_v2_m1.exons.binom_tests.txt",quote=FALSE,row.names=FALSE,sep="\t");

write.table(data3_tests,file="/Users/kraigrs/Wittkopp/Graze/SRR389095.mate1.berlin-updated-exonic-regions.bowtie_v3_m1.exons.binom_tests.txt",quote=FALSE,row.names=FALSE,sep="\t");

##########

par(mfrow=c(2,2));

boxplot(log2(berlin_ref_allele/c1674_alt_allele) ~ neighbors, data = subset(data0_complete,is.finite(log2(berlin_ref_allele/c1674_alt_allele))),varwidth = TRUE,xlab="",ylab="",main="multiple references",xlim=c(0,34),ylim=c(-8,13),col=rgb(0,0,0,0.3),pars = list(outpch=19,cex=0.2));

boxplot(log2(ref_allele/alt_allele) ~ neighbors, data = subset(data1_complete,is.finite(log2(ref_allele/alt_allele))),varwidth = TRUE,xlab="",ylab="",main="single reference (1mm)",xlim=c(0,34),ylim=c(-8,13),col=rgb(0,0,0,0.3),pars = list(outpch=19,cex=0.2));

boxplot(log2(ref_allele/alt_allele) ~ neighbors, data = subset(data2_complete,is.finite(log2(ref_allele/alt_allele))),varwidth = TRUE,xlab="",ylab="",main="single reference (2mm)",xlim=c(0,34),ylim=c(-8,13),col=rgb(0,0,0,0.3),pars = list(outpch=19,cex=0.2));

boxplot(log2(ref_allele/alt_allele) ~ neighbors, data = subset(data3_complete,is.finite(log2(ref_allele/alt_allele))),varwidth = TRUE,xlab="",ylab="",main="single reference (3mm)",xlim=c(0,34),ylim=c(-8,13),col=rgb(0,0,0,0.3),pars = list(outpch=19,cex=0.2));

summary(log2(data0_complete$berlin_ref_allele/data0_complete$c1674_alt_allele));
summary(log2(data1_complete$ref_allele/data1_complete$alt_allele));
summary(log2(data2_complete$ref_allele/data2_complete$alt_allele));
summary(log2(data3_complete$ref_allele/data3_complete$alt_allele));

par(mfrow=c(2,2));

boxplot( ref_allele/(ref_allele+alt_allele) ~ neighbors_plus, data = data1_complete,varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,29),ylim=c(0,1),pars = list(outpch=19,cex=0.2));
abline(h=0.5,lty=2,col="red");

boxplot(ref_allele/(ref_allele+alt_allele) ~ neighbors_plus, data = data2_complete,varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,29),ylim=c(0,1),pars = list(outpch=19,cex=0.2));
abline(h=0.5,lty=2,col="red");

boxplot(ref_allele/(ref_allele+alt_allele) ~ neighbors_plus, data = data3_complete,varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,29),ylim=c(0,1),pars = list(outpch=19,cex=0.2));
abline(h=0.5,lty=2,col="red");

boxplot(berlin_ref_allele/(berlin_ref_allele+c1674_alt_allele) ~ neighbors_plus, data = data0_complete,varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,29),ylim=c(0,1),pars = list(outpch=19,cex=0.2));
abline(h=0.5,lty=2,col="red");

props <- c(65364/220928,32402/291700,22768/339658,10038/327693);
sig <- c(65364,32402,22768,10038);
tot <- c(220928,291700,339658,327693);

#barplot(props,names.arg=c(1,2,3,0),ylim=c(0,1),xlab="Number of mismatches",ylab="Proportion of AI",main="Comparison of single and multiple genomes");

pie(c(sig[1],tot[1]-sig[1]),col=c("gray","white"),labels="");
pie(c(sig[2],tot[2]-sig[2]),col=c("gray","white"),labels="");
pie(c(sig[3],tot[3]-sig[3]),col=c("gray","white"),labels="");
pie(c(sig[4],tot[4]-sig[4]),col=c("gray","white"),labels="");


barplot(props,names.arg=c(1,2,3,0),ylim=c(0,1),xlab="",ylab="",main="");

barplot(c(0.98,0.08),names.arg=c("Simulated","Real"),ylim=c(0,1),xlab="",ylab="",main="");


###### mappability ######

data1_tests <- read.table("/Users/kraigrs/Wittkopp/Graze/SRR389095.mate1.berlin-updated-exonic-regions.bowtie_v1_m1.SNPs.binom_tests.txt",header=TRUE,sep="\t");

data2_tests <- read.table("/Users/kraigrs/Wittkopp/Graze/SRR389095.mate1.berlin-updated-exonic-regions.bowtie_v2_m1.SNPs.binom_tests.txt",header=TRUE,sep="\t");

data3_tests <- read.table("/Users/kraigrs/Wittkopp/Graze/SRR389095.mate1.berlin-updated-exonic-regions.bowtie_v3_m1.SNPs.binom_tests.txt",header=TRUE,sep="\t");

data0_tests <- read.table("/Users/kraigrs/Wittkopp/Graze/SRR389095.mate1.berlin_c1674_fsa_masked.bowtie_v0_m1.SNPs.binom_tests.txt",header=TRUE,sep="\t");

neighbors_plus <- data0_tests$neighbors+1;
data0_tests <- cbind(data0_tests,neighbors_plus);

neighbors_plus <- data1_tests$neighbors+1;
data1_tests <- cbind(data1_tests,neighbors_plus);

neighbors_plus <- data2_tests$neighbors+1;
data2_tests <- cbind(data2_tests,neighbors_plus);

neighbors_plus <- data3_tests$neighbors+1;
data3_tests <- cbind(data3_tests,neighbors_plus);


map1 <- read.table("/Users/kraigrs/Wittkopp/Graze/berlin-updated-exonic-regions.SNPs.l36_m1.mappability.txt",header=TRUE,sep="\t");
map2 <- read.table("/Users/kraigrs/Wittkopp/Graze/berlin-updated-exonic-regions.SNPs.l36_m2.mappability.txt",header=TRUE,sep="\t");
map3 <- read.table("/Users/kraigrs/Wittkopp/Graze/berlin-updated-exonic-regions.SNPs.l36_m3.mappability.txt",header=TRUE,sep="\t");

berlin_map <- read.table("/Users/kraigrs/Wittkopp/Graze/berlin-updated-exonic-regions_fsa_masked.SNPs.mappability.txt",header=TRUE,sep="\t");
c1674_map <- read.table("/Users/kraigrs/Wittkopp/Graze/c1674-updated-exonic-regions_fsa_masked.SNPs.mappability.txt",header=TRUE,sep="\t");
map0 <- merge(berlin_map,c1674_map,by.x=c("locus","position"),by.y=c("locus","position"));

#summary(data0_tests_map0$neighbors);
#summary(data1_tests_map1$neighbors);
#summary(data2_tests_map2$neighbors);
#summary(data3_tests_map3$neighbors);


# 1 mismatch

data1_tests_map1 <- merge(data1_tests,map1,by.x=c("chr","pos"),by.y=c("locus","position"));

write.table(data1_tests_map1,file="/Users/kraigrs/Desktop/exons/real_SNPs_1mm.txt",quote=F,sep="\t",row.names=F,col.names=F);

temp <- subset(data1_tests_map1,neighbors<1);

# proportion of differentiating sites with perfect mappability and equal allelic abundance
nrow(subset(temp,qvals >= 0.05 & sum/length == 1))/nrow(subset(temp,sum/length == 1))*100;  

nrow(subset(temp,qvals < 0.05 & sum/length != 1))/nrow(subset(temp,qvals < 0.05))*100; # % unequal ASRA and imperfect mappability 
nrow(subset(temp,qvals >= 0.05 & sum/length != 1))/nrow(subset(temp,qvals >= 0.05))*100; # % equal ASRA and imperfect mappability

perfect <- subset(temp,sum/length == 1);
imperfect <- subset(temp,sum/length != 1);

perfect_hist <- hist(perfect$ref_allele/(perfect$ref_allele+perfect$alt_allele),breaks=20);
imperfect_hist <- hist(imperfect$ref_allele/(imperfect$ref_allele+imperfect$alt_allele),breaks=20);
mat <- cbind(perfect_hist$counts/nrow(perfect),imperfect_hist$counts/nrow(imperfect));
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),xlab="",ylab="",main="",ylim=c(0,0.3),col=c("white","grey"));

nonsig <- subset(temp,qvals >= 0.05);
sig <- subset(temp,qvals < 0.05);
nrow(subset(sig,sum/length != 1 ))/nrow(sig)*100;

nonsig_hist <- hist( (nonsig$sum/nonsig$length) ,breaks=seq(0,1,0.05));
sig_hist <- hist( (sig$sum/sig$length) ,breaks=seq(0,1,0.05));
mat <- cbind(nonsig_hist$counts/nrow(nonsig),sig_hist$counts/nrow(sig));
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),xlab="",ylab="",main="",col=c("white","grey"));

perfect <- subset(data1_tests_map1,sum/length == 1);
imperfect <- subset(data1_tests_map1,sum/length != 1);

perfect_hist <- hist(perfect$ref_allele/(perfect$ref_allele+perfect$alt_allele),breaks=seq(0,1,0.05));
imperfect_hist <- hist(imperfect$ref_allele/(imperfect$ref_allele+imperfect$alt_allele),breaks=seq(0,1,0.05));
mat <- cbind(perfect_hist$counts/nrow(perfect),imperfect_hist$counts/nrow(imperfect));

pdf(file="/Users/kraigrs/Desktop/plots/real_1mm_barplot_by_mapp.pdf");
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),ylim=c(0,1),xlab="",ylab="",main="",col=c("white","grey"));
#par(new=TRUE);
#points(seq(0.05,1,0.05),mat[,1],ylim=c(0,1),col="black",type="b");
#par(new=TRUE);
#points(seq(0.05,1,0.05),mat[,2],ylim=c(0,1),col="grey",type="b");
dev.off();

pdf(file="/Users/kraigrs/Desktop/plots/real_1mm_pie_by_mapp.pdf");
pie(nrow(perfect),nrow(imperfect),labels="");
dev.off();

pdf(file="/Users/kraigrs/Desktop/plots/real_1mm_boxplot_perfect.pdf");
boxplot(ref_allele/(ref_allele+alt_allele) ~ neighbor_plus,data = perfect, varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,29),ylim=c(0,1),pars = list(outpch=19,cex=0.4,col=rgb(0,0,0,0.4)));
abline(h = 0.5,col="red",lty="dashed");
dev.off();

pdf(file="/Users/kraigrs/Desktop/plots/real_1mm_boxplot_imperfect.pdf");
boxplot(ref_allele/(ref_allele+alt_allele) ~ neighbor_plus,data = imperfect, varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,29),ylim=c(0,1),pars = list(outpch=19,cex=0.4,col=rgb(0,0,0,0.4)));
abline(h = 0.5,col="red",lty="dashed");
dev.off();


d1_perfect <- density(perfect$ref_allele/(perfect$ref_allele+perfect$alt_allele));
d1_imperfect <- density(imperfect$ref_allele/(imperfect$ref_allele+imperfect$alt_allele));

plot(d1_perfect,xlim=c(0,1),ylim=c(0,10),main="",col="blue",xlab="");
par(new=TRUE);
plot(d1_imperfect,xlim=c(0,1),ylim=c(0,10),main="",col="red",xlab="");

# 2 mismatches

data2_tests_map2 <- merge(data2_tests,map2,by.x=c("chr","pos"),by.y=c("locus","position"));

write.table(data2_tests_map2,file="/Users/kraigrs/Desktop/exons/real_SNPs_2mm.txt",quote=F,sep="\t",row.names=F,col.names=F);

temp <- subset(data2_tests_map2,neighbors<2);

# proportion of differentiating sites with perfect mappability and equal allelic abundance
nrow(subset(temp,qvals >= 0.05 & sum/length == 1))/nrow(subset(temp,sum/length == 1))*100;  

nrow(subset(temp,qvals < 0.05 & sum/length != 1))/nrow(subset(temp,qvals < 0.05))*100; # % unequal ASRA and imperfect mappability 
nrow(subset(temp,qvals >= 0.05 & sum/length != 1))/nrow(subset(temp,qvals >= 0.05))*100; # % equal ASRA and imperfect mappability

perfect <- subset(temp,sum/length == 1);
imperfect <- subset(temp,sum/length != 1);

perfect_hist <- hist(perfect$ref_allele/(perfect$ref_allele+perfect$alt_allele),breaks=20);
imperfect_hist <- hist(imperfect$ref_allele/(imperfect$ref_allele+imperfect$alt_allele),breaks=20);
mat <- cbind(perfect_hist$counts/nrow(perfect),imperfect_hist$counts/nrow(imperfect));
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),xlab="",ylab="",main="",ylim=c(0,0.3),col=c("white","grey"));

nonsig <- subset(temp,qvals >= 0.05);
sig <- subset(temp,qvals < 0.05);
nrow(subset(sig,sum/length != 1 ))/nrow(sig)*100;

nonsig_hist <- hist( (nonsig$sum/nonsig$length) ,breaks=seq(0,1,0.05));
sig_hist <- hist( (sig$sum/sig$length) ,breaks=seq(0,1,0.05));
mat <- cbind(nonsig_hist$counts/nrow(nonsig),sig_hist$counts/nrow(sig));
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),xlab="",ylab="",main="",col=c("white","grey"));

perfect <- subset(data2_tests_map2,sum/length == 1);
imperfect <- subset(data2_tests_map2,sum/length != 1);
boxplot( ref_allele/(ref_allele+alt_allele) ~ neighbors_plus, data = perfect,varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,29),ylim=c(0,1),pars = list(outpch=19,cex=0.2));
abline(h=0.5,lty=2,col="red");

d2_perfect <- density(perfect$ref_allele/(perfect$ref_allele+perfect$alt_allele));
d2_imperfect <- density(imperfect$ref_allele/(imperfect$ref_allele+imperfect$alt_allele));

plot(d2_perfect,xlim=c(0,1),ylim=c(0,10),main="",col="blue",xlab="");
par(new=TRUE);
plot(d2_imperfect,xlim=c(0,1),ylim=c(0,10),main="",col="red",xlab="");

# 3 mismatches

data3_tests_map3 <- merge(data3_tests,map3,by.x=c("chr","pos"),by.y=c("locus","position"));

write.table(data3_tests_map3,file="/Users/kraigrs/Desktop/exons/real_SNPs_3mm.txt",quote=F,sep="\t",row.names=F,col.names=F);

temp <- subset(data3_tests_map3,neighbors<3);

# proportion of differentiating sites with perfect mappability and equal allelic abundance
nrow(subset(temp,qvals >= 0.05 & sum/length == 1))/nrow(subset(temp,sum/length == 1))*100;

nrow(subset(temp,qvals < 0.05 & sum/length != 1))/nrow(subset(temp,qvals < 0.05))*100; # % unequal ASRA and imperfect mappability 
nrow(subset(temp,qvals >= 0.05 & sum/length != 1))/nrow(subset(temp,qvals >= 0.05))*100; # % equal ASRA and imperfect mappability

perfect <- subset(temp,sum/length == 1);
imperfect <- subset(temp,sum/length != 1);

perfect_hist <- hist(perfect$ref_allele/(perfect$ref_allele+perfect$alt_allele),breaks=20);
imperfect_hist <- hist(imperfect$ref_allele/(imperfect$ref_allele+imperfect$alt_allele),breaks=20);
mat <- cbind(perfect_hist$counts/nrow(perfect),imperfect_hist$counts/nrow(imperfect));
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),xlab="",ylab="",main="",ylim=c(0,0.3),col=c("white","grey"));

nonsig <- subset(temp,qvals >= 0.05);
sig <- subset(temp,qvals < 0.05);
nrow(subset(sig,sum/length != 1 ))/nrow(sig)*100;

nonsig_hist <- hist( (nonsig$sum/nonsig$length) ,breaks=seq(0,1,0.05));
sig_hist <- hist( (sig$sum/sig$length) ,breaks=seq(0,1,0.05));
mat <- cbind(nonsig_hist$counts/nrow(nonsig),sig_hist$counts/nrow(sig));
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),xlab="",ylab="",main="",col=c("white","grey"));

perfect <- subset(data3_tests_map3,sum/length == 1);
imperfect <- subset(data3_tests_map3,sum/length != 1);
boxplot( ref_allele/(ref_allele+alt_allele) ~ neighbors_plus, data = perfect,varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,29),ylim=c(0,1),pars = list(outpch=19,cex=0.2));
abline(h=0.5,lty=2,col="red");

d3_perfect <- density(perfect$ref_allele/(perfect$ref_allele+perfect$alt_allele));
d3_imperfect <- density(imperfect$ref_allele/(imperfect$ref_allele+imperfect$alt_allele));

plot(d3_perfect,xlim=c(0,1),ylim=c(0,10),main="",col="blue",xlab="");
par(new=TRUE);
plot(d3_imperfect,xlim=c(0,1),ylim=c(0,10),main="",col="red",xlab="");

# 0 mismatches

data0_tests_map0 <- merge(data0_tests,map0,by.x=c("chr","pos"),by.y=c("locus","position"));

write.table(data0_tests_map0,file="/Users/kraigrs/Desktop/exons/real_SNPs_0mm.txt",quote=F,sep="\t",row.names=F,col.names=F);

temp <- data0_tests_map0;

# proportion of differentiating sites with perfect mappability and equal allelic abundance
nrow(subset(temp,qvals >= 0.05 & ((sum.x/length.x)+(sum.y/length.y)) == 2))/nrow(subset(temp,((sum.x/length.x)+(sum.y/length.y)) == 2))*100;

nrow(subset(temp,qvals < 0.05 & ((sum.x/length.x)+(sum.y/length.y)) < 2))/nrow(subset(temp,qvals < 0.05))*100; # % unequal ASRA and imperfect mappability 
nrow(subset(temp,qvals >= 0.05 & ((sum.x/length.x)+(sum.y/length.y)) < 2))/nrow(subset(temp,qvals >= 0.05))*100; # % equal ASRA and imperfect mappability

perfect <- subset(temp,((sum.x/length.x)+(sum.y/length.y)) == 2);
imperfect <- subset(temp,((sum.x/length.x)+(sum.y/length.y)) < 2);

perfect_hist <- hist(perfect$berlin_ref_allele/(perfect$berlin_ref_allele+perfect$c1674_alt_allele),breaks=20);
imperfect_hist <- hist(imperfect$berlin_ref_allele/(imperfect$berlin_ref_allele+imperfect$c1674_alt_allele),breaks=20);
mat <- cbind(perfect_hist$counts/nrow(perfect),imperfect_hist$counts/nrow(imperfect));
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),ylim=c(0,0.3),xlab="",ylab="",main="",col=c("white","grey"));

nonsig <- subset(temp,qvals >= 0.05);
sig <- subset(temp,qvals < 0.05);
nrow(subset(sig,sum.x/length.x != 1 | sum.y/length.y != 1))/nrow(sig)*100;

nonsig_hist <- hist( (nonsig$sum.x/nonsig$length.x)+(nonsig$sum.y/nonsig$length.y) ,breaks=seq(0,2,0.1));
sig_hist <- hist( (sig$sum.x/sig$length.x)+(sig$sum.y/sig$length.y) ,breaks=seq(0,2,0.1));
mat <- cbind(nonsig_hist$counts/nrow(nonsig),sig_hist$counts/nrow(sig));
barplot(t(mat),beside=TRUE,names.arg=seq(0.1,2,0.1),xlab="",ylab="",main="",col=c("white","grey"));

perfect <- subset(data0_tests_map0,((sum.x/length.x)+(sum.y/length.y)) == 2);
imperfect <- subset(data0_tests_map0,((sum.x/length.x)+(sum.y/length.y)) != 2);
boxplot( berlin_ref_allele/(berlin_ref_allele+c1674_alt_allele) ~ neighbors_plus, data = perfect,varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,29),ylim=c(0,1),pars = list(outpch=19,cex=0.2));
abline(h=0.5,lty=2,col="red");

d0_perfect <- density(perfect$berlin_ref_allele/(perfect$berlin_ref_allele+perfect$c1674_alt_allele));
d0_imperfect <- density(imperfect$berlin_ref_allele/(imperfect$berlin_ref_allele+imperfect$c1674_alt_allele));

plot(d0_perfect,xlim=c(0,1),ylim=c(0,10),main="",col="blue",xlab="");
par(new=TRUE);
plot(d0_imperfect,xlim=c(0,1),ylim=c(0,10),main="",col="red",xlab="");

######### make plots comparing approaches ############

key <- read.table("/Users/kraigrs/Wittkopp/Graze/key",header=TRUE,sep="\t");

perfect0 <- subset(data0_tests_map0,((sum.x/length.x)+(sum.y/length.y)) == 2);

temp1 <- subset(data1_tests_map1,neighbors<1);
perfect1 <- subset(temp1,sum/length == 1);

temp2 <- subset(data2_tests_map2,neighbors<2);
perfect2 <- subset(temp2,sum/length == 1);

temp3 <- subset(data3_tests_map3,neighbors<3);
#temp3 <- data3_tests_map3;
perfect3 <- subset(temp3,sum/length == 1);


x <- merge(perfect0,key,by.x=c("chr","pos"),by.y=c("locus","fsa"));
y <- merge(perfect1,key,by.x=c("chr","pos"),by.y=c("locus","mel"));
compare <- merge(x,y,by.x=c("chr","mel"),by.y=c("chr","pos"));
nrow(compare);
s_ns <- subset(compare,qvals.x < 0.05 & qvals.y >= 0.05);
ns_s <- subset(compare,qvals.x >= 0.05 & qvals.y < 0.05);
tot <- nrow(s_ns)+nrow(ns_s);
a <- nrow(subset(s_ns,ref_allele/(ref_allele+alt_allele) > berlin_ref_allele/(berlin_ref_allele+c1674_alt_allele)));
b <- nrow(subset(ns_s,ref_allele/(ref_allele+alt_allele) > berlin_ref_allele/(berlin_ref_allele+c1674_alt_allele)));
(a+b)/tot;

pdf(file="/Users/kraigrs/Desktop/new_plots/qqplot_real_data_qvals_3mm_unfiltered.pdf");
qqplot(-log10(compare$qvals.x),-log10(compare$qvals.y),xlab="-log10(qvals) parental genomes",ylab="-log10(qvals) single genome",main="");
abline(a=0,b=1,col="red",lty=2);
dev.off();

pdf(file="/Users/kraigrs/Desktop/new_plots/compare_real_data_qvals_3mm_unfiltered.pdf");

plot( compare$berlin_ref_allele/(compare$berlin_ref_allele+compare$c1674_alt_allele) , compare$ref_allele/(compare$ref_allele+compare$alt_allele) ,xlim=c(0,1),ylim=c(0,1),type="n",pch=19,cex=0.5,col=rgb(0,0,0,0.05),xlab="Maternal and paternal genomes, 0 mismatches",ylab="Single genome, 3 mismatches (unfiltered)",main="Compare proportion of reference allele\nfor perfect mappability");

ns_ns <- subset(compare,qvals.x >= 0.05 & qvals.y >= 0.05);
nrow(ns_ns);
points(ns_ns$berlin_ref_allele/(ns_ns$berlin_ref_allele+ns_ns$c1674_alt_allele),ns_ns$ref_allele/(ns_ns$ref_allele+ns_ns$alt_allele),pch=19,cex=0.5,col=rgb(0,0,0,0.05));

s_ns <- subset(compare,qvals.x < 0.05 & qvals.y >= 0.05);
nrow(s_ns);
points(s_ns$berlin_ref_allele/(s_ns$berlin_ref_allele+s_ns$c1674_alt_allele),s_ns$ref_allele/(s_ns$ref_allele+s_ns$alt_allele),pch=19,cex=0.5,col=rgb(1,0,0,0.6));

ns_s <- subset(compare,qvals.x >= 0.05 & qvals.y < 0.05);
nrow(ns_s);
points(ns_s$berlin_ref_allele/(ns_s$berlin_ref_allele+ns_s$c1674_alt_allele),ns_s$ref_allele/(ns_s$ref_allele+ns_s$alt_allele),pch=19,cex=0.5,col=rgb(0,0,1,0.6));

s_s <- subset(compare,qvals.x < 0.05 & qvals.y < 0.05);
nrow(s_s);
points(s_s$berlin_ref_allele/(s_s$berlin_ref_allele+s_s$c1674_alt_allele),s_s$ref_allele/(s_s$ref_allele+s_s$alt_allele),pch=19,cex=0.5,col=rgb(1,0,1,0.6));

dev.off();

ns_ns <- subset(compare,pvals.x >= 0.05 & pvals.y >= 0.05);
points(ns_ns$berlin_ref_allele/(ns_ns$berlin_ref_allele+ns_ns$c1674_alt_allele),ns_ns$ref_allele/(ns_ns$ref_allele+ns_ns$alt_allele),pch=19,cex=0.5,col=rgb(0,0,0,0.2));

s_ns <- subset(compare,pvals.x < 0.05 & pvals.y >= 0.05);
points(s_ns$berlin_ref_allele/(s_ns$berlin_ref_allele+s_ns$c1674_alt_allele),s_ns$ref_allele/(s_ns$ref_allele+s_ns$alt_allele),pch=19,cex=0.5,col=rgb(1,0,0,0.2));

ns_s <- subset(compare,pvals.x >= 0.05 & pvals.y < 0.05);
points(ns_s$berlin_ref_allele/(ns_s$berlin_ref_allele+ns_s$c1674_alt_allele),ns_s$ref_allele/(ns_s$ref_allele+ns_s$alt_allele),pch=19,cex=0.5,col=rgb(0,0,1,0.2));

s_s <- subset(compare,pvals.x < 0.05 & pvals.y < 0.05);
points(s_s$berlin_ref_allele/(s_s$berlin_ref_allele+s_s$c1674_alt_allele),s_s$ref_allele/(s_s$ref_allele+s_s$alt_allele),pch=19,cex=0.5,col=rgb(1,0,1,0.2));

###########
data <- merge(data0_tests,map,by.x=c("chr","pos"),by.y=c("locus","position"));
#data <- merge(data0_complete,map,by.x=c("chr","pos"),by.y=c("locus","position"));

test <- subset(data,log2(berlin_ref_allele/c1674_alt_allele) > 7.5);
compare <- merge(test,data3_complete,by.x=c("chr","pos"),by.y=c("chr","pos"));

par(mfrow=c(1,2));

hist(log2(compare$ref_allele/compare$alt_allele),main="Single reference (3mm)",xlab="log2(ASE)");
hist(log2(compare$berlin_ref_allele/compare$c1674_alt_allele),main="Multiple references (0mm)",xlab="log2(ASE)");

par(mfrow=c(1,2));

plot(compare$berlin_ref_allele,compare$ref_allele,ylab="Reference allele from single reference (3mm)",xlab="Reference allele from multiple references (0mm)",pch=19,col=rgb(0,0,0,0.3),cex=0.5);
plot(compare$c1674_alt_allele,compare$alt_allele,ylab="Alternative allele from single reference (3mm)",xlab="Alternative allele from multiple references (0mm)",pch=19,col=rgb(0,0,0,0.3),cex=0.5);

plot(log2(compare$berlin_ref_allele/compare$c1674_alt_allele),log2(compare$ref_allele/compare$alt_allele),xlim=c(2.5,13),ylim=c(2.5,13),xlab="log2(ASE) multiple references (0mm)",ylab="log2(ASE) single reference (3mm)");
abline(a=0,b=1,col="red",lty=2);

plot(log2(data$berlin_ref_allele/data$c1674_alt_allele),log2((data$sum.x/data$length.x)/(data$sum.y/data$length.y)),pch=19,cex=0.3,col=rgb(0,0,0,0.3));

check <- subset(data,
	abs(log2(berlin_ref_allele/c1674_alt_allele)) > log2(19)
	);

perfect <- subset(check,
	log2((sum.x/length.x)/(sum.y/length.y)) == 0
	);
	
imperfect <- subset(check,
	log2((sum.x/length.x)/(sum.y/length.y)) != 0
	);

summary(log2((check$sum.x/check$length.x)/(check$sum.y/check$length.y)));

hist(log2(check$berlin_ref_allele/check$c1674_alt_allele));
hist(check$berlin_ref_allele/check$c1674_alt_allele);

nrow(subset(check,log2((sum.x/length.x)/(sum.y/length.y))==0 ))/nrow(check);

###########################################################################
# explore some properties of density of variable sites and read abundance #
###########################################################################

quantile(data0_complete$neighbors,probs=0.99);

extreme <- subset(data0_complete,neighbors>11);
real <- subset(data0_complete,neighbors<=11);

par(mfrow=c(1,2));
hist(extreme$berlin_ref_allele/(extreme$berlin_ref_allele+extreme$c1674_alt_allele),breaks=20,freq=FALSE,xlab="mel/(mel+sim)",main="Variable site density > 11");
hist(real$berlin_ref_allele/(real$berlin_ref_allele+real$c1674_alt_allele),breaks=20,freq=FALSE,xlab="mel/(mel+sim)",main="Variable site density <= 11");

par(mfrow=c(1,2));
hist(log2(extreme$berlin_ref_allele+extreme$c1674_alt_allele),breaks=20,freq=FALSE,xlab="log2(mel+sim)",main="Variable site density > 11",xlim=c(4,8));
hist(log2(real$berlin_ref_allele+real$c1674_alt_allele),breaks=40,freq=FALSE,xlab="log2(mel+sim)",main="Variable site density <= 11",xlim=c(4,8));




