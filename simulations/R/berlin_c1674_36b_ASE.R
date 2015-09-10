data1 <- read.table("/Users/kraigrs/Wittkopp/Graze/simulation/berlin_c1674.tiled_36bp.berlin-updated-exonic-regions.bowtie_v1_m1.SNPs.txt",header=TRUE,sep="\t");

data2 <- read.table("/Users/kraigrs/Wittkopp/Graze/simulation/berlin_c1674.tiled_36bp.berlin-updated-exonic-regions.bowtie_v2_m1.SNPs.txt",header=TRUE,sep="\t");

data3 <- read.table("/Users/kraigrs/Wittkopp/Graze/simulation/berlin_c1674.tiled_36bp.berlin-updated-exonic-regions.bowtie_v3_m1.SNPs.txt",header=TRUE,sep="\t");

data0_ref <- read.table("/Users/kraigrs/Wittkopp/Graze/simulation/berlin_c1674.tiled_36bp.berlin-updated-exonic-regions_fsa_masked.bowtie_v0_m1.SNPs.txt",header=TRUE,sep="\t");

data0_alt <- read.table("/Users/kraigrs/Wittkopp/Graze/simulation/berlin_c1674.tiled_36bp.c1674-updated-exonic-regions_fsa_masked.bowtie_v0_m1.SNPs.txt",header=TRUE,sep="\t");

data0 <- merge(data0_ref,data0_alt,by.x=c("chr","pos"),by.y=c("chr","pos"));




summary(data0$neighbors.x);
summary(data1$neighbors);
summary(data2$neighbors);
summary(data3$neighbors);

par(mfrow=c(2,2));

hist(data1$neighbors,breaks=100);
hist(data2$neighbors,breaks=100);
hist(data3$neighbors,breaks=100);
hist(data0$neighbors.x,breaks=100);

quantile(data1$neighbors,probs=c(0.95,0.975,0.999));
quantile(data2$neighbors,probs=c(0.95,0.975,0.999));
quantile(data3$neighbors,probs=c(0.95,0.975,0.999));
quantile(data0$neighbors.x,probs=c(0.95,0.975,0.999));

d0 <- density(data0_complete$berlin_ref_allele/(data0_complete$berlin_ref_allele+data0_complete$c1674_alt_allele));
d1 <- density(data1_complete$ref_allele/(data1_complete$ref_allele+data1_complete$alt_allele));
d2 <- density(data2_complete$ref_allele/(data2_complete$ref_allele+data2_complete$alt_allele));
d3 <- density(data3_complete$ref_allele/(data3_complete$ref_allele+data3_complete$alt_allele));

par(mfrow=c(2,2));

plot(d1,xlim=c(0,1),ylim=c(0,8),main="single reference (1mm)");
plot(d2,xlim=c(0,1),ylim=c(0,8),main="single reference (2mm)");
plot(d3,xlim=c(0,1),ylim=c(0,8),main="single reference (3mm)");
plot(d0,xlim=c(0,1),ylim=c(0,8),main="multiple references");

nrow(subset(data0_complete,log2(berlin_ref_allele/c1674_alt_allele) != 0))/nrow(data0_complete);
nrow(subset(data3_complete,log2(ref_allele/alt_allele) != 0))/nrow(data3_complete);
nrow(subset(data2_complete,log2(ref_allele/alt_allele) != 0))/nrow(data2_complete);
nrow(subset(data1_complete,log2(ref_allele/alt_allele) != 0))/nrow(data1_complete);


##########

neighbors_plus <- data1$neighbors+1;
data1 <- cbind(data1,neighbors_plus);

neighbors_plus <- data2$neighbors+1;
data2 <- cbind(data2,neighbors_plus);

neighbors_plus <- data3$neighbors+1;
data3 <- cbind(data3,neighbors_plus);

neighbors_plus <- data0$neighbors.x+1;
data0 <- cbind(data0,neighbors_plus);

#par(mfrow=c(2,2));

boxplot( ref_allele/(ref_allele+alt_allele) ~ neighbors_plus, data = data1,varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,15),ylim=c(0,1),pars = list(outpch=19,cex=0.2));
abline(h=0.5,lty=2,col="red");

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_36b_v1_box.pdf");
boxplot( ref_allele/(ref_allele+alt_allele) ~ neighbors_plus, data = data1,varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,15),ylim=c(0,1),pars = list(outpch=19,outcex=0.2));
abline(h=0.5,lty=2,col="red");
dev.off();


postscript(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_36b_v1_box.ps",paper="special",width=7,height=7);
boxplot( ref_allele/(ref_allele+alt_allele) ~ neighbors_plus, data = data1,varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,15),ylim=c(0,1),pars = list(outpch=19,outcex=0.2));
abline(h=0.5,lty=2,col="red");
dev.off();

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_36b_v2_box.pdf");
boxplot(ref_allele/(ref_allele+alt_allele) ~ neighbors_plus, data = data2,varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,15),ylim=c(0,1),pars = list(outpch=19,outcex=0.2));
abline(h=0.5,lty=2,col="red");
dev.off();

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_36b_v3_box.pdf");
boxplot(ref_allele/(ref_allele+alt_allele) ~ neighbors_plus, data = data3,varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,15),ylim=c(0,1),pars = list(outpch=19,outcex=0.2));
abline(h=0.5,lty=2,col="red");
dev.off();

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_36b_v0_box.pdf");
boxplot(ref_allele.x/(ref_allele.x+ref_allele.y) ~ neighbors_plus, data = data0,varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,15),ylim=c(0,1),pars = list(outpch=19,outcex=0.2));
abline(h=0.5,lty=2,col="red");
dev.off();


pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_36b_v1_pie.pdf");
pie(c(nrow(subset(data1,ref_allele/(ref_allele+alt_allele) != 0.5)),nrow(subset(data1,ref_allele/(ref_allele+alt_allele) == 0.5))),col=c("gray","white"),labels="");
dev.off();

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_36b_v2_pie.pdf");
pie(c(nrow(subset(data2,ref_allele/(ref_allele+alt_allele) != 0.5)),nrow(subset(data2,ref_allele/(ref_allele+alt_allele) == 0.5))),col=c("gray","white"),labels="");
dev.off();

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_36b_v3_pie.pdf");
pie(c(nrow(subset(data3,ref_allele/(ref_allele+alt_allele) != 0.5)),nrow(subset(data3,ref_allele/(ref_allele+alt_allele) == 0.5))),col=c("gray","white"),labels="");
dev.off();

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_36b_v0_pie.pdf");
pie(c(nrow(subset(data0,ref_allele.x/(ref_allele.x+ref_allele.y) != 0.5)),nrow(subset(data0,ref_allele.x/(ref_allele.x+ref_allele.y) == 0.5))),col=c("gray","white"),labels="");
dev.off();

###### mappability ######


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

data1_map1 <- merge(data1,map1,by.x=c("chr","pos"),by.y=c("locus","position"));
write.table(data1_map1,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_SNPs_36b_1mm.txt",quote=F,sep="\t",row.names=F,col.names=F);
unbiased <- subset(data1_map1,neighbors<1);
write.table(unbiased,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_SNPs_36b_1mm_unbiased.txt",quote=F,sep="\t",row.names=F,col.names=F);
perfect <- subset(unbiased,sum/length == 1);
write.table(perfect,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_SNPs_36b_1mm_perfect.txt",quote=F,sep="\t",row.names=F,col.names=F);
nrow(unbiased);
nrow(perfect);

nrow(subset(data1_map1,ref_allele/(ref_allele+alt_allele) != 0.5 & sum/length != 1))/nrow(subset(data1_map1,ref_allele/(ref_allele+alt_allele) != 0.5));

perfect <- subset(data1_map1,sum/length == 1);
imperfect <- subset(data1_map1,sum/length != 1);
perfect_hist <- hist(perfect$ref_allele/(perfect$ref_allele+perfect$alt_allele),breaks=seq(0,1,0.05));
imperfect_hist <- hist(imperfect$ref_allele/(imperfect$ref_allele+imperfect$alt_allele),breaks=seq(0,1,0.05));
mat <- cbind(perfect_hist$counts/nrow(perfect),imperfect_hist$counts/nrow(imperfect));

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_1mm_barplot_by_mapp.pdf");
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),ylim=c(0,1),xlab="",ylab="",main="",col=c("white","grey"));
dev.off();


data2_map2 <- merge(data2,map2,by.x=c("chr","pos"),by.y=c("locus","position"));
write.table(data2_map2,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_SNPs_36b_2mm.txt",quote=F,sep="\t",row.names=F,col.names=F);
unbiased <- subset(data2_map2,neighbors<2);
write.table(unbiased,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_SNPs_36b_2mm_unbiased.txt",quote=F,sep="\t",row.names=F,col.names=F);
perfect <- subset(unbiased,sum/length == 1);
write.table(perfect,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_SNPs_36b_2mm_perfect.txt",quote=F,sep="\t",row.names=F,col.names=F);
nrow(unbiased);
nrow(perfect);

nrow(subset(data2_map2,ref_allele/(ref_allele+alt_allele) != 0.5 & sum/length != 1))/nrow(subset(data2_map2,ref_allele/(ref_allele+alt_allele) != 0.5));

perfect <- subset(data2_map2,sum/length == 1);
imperfect <- subset(data2_map2,sum/length != 1);
perfect_hist <- hist(perfect$ref_allele/(perfect$ref_allele+perfect$alt_allele),breaks=seq(0,1,0.05));
imperfect_hist <- hist(imperfect$ref_allele/(imperfect$ref_allele+imperfect$alt_allele),breaks=seq(0,1,0.05));
mat <- cbind(perfect_hist$counts/nrow(perfect),imperfect_hist$counts/nrow(imperfect));

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_2mm_barplot_by_mapp.pdf");
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),ylim=c(0,1),xlab="",ylab="",main="",col=c("white","grey"));
dev.off();


data3_map3 <- merge(data3,map3,by.x=c("chr","pos"),by.y=c("locus","position"));
write.table(data3_map3,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_SNPs_36b_3mm.txt",quote=F,sep="\t",row.names=F,col.names=F);
unbiased <- subset(data3_map3,neighbors<3);
write.table(unbiased,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_SNPs_36b_3mm_unbiased.txt",quote=F,sep="\t",row.names=F,col.names=F);
perfect <- subset(unbiased,sum/length == 1);
write.table(perfect,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_SNPs_36b_3mm_perfect.txt",quote=F,sep="\t",row.names=F,col.names=F);
nrow(unbiased);
nrow(perfect);

nrow(subset(data3_map3,ref_allele/(ref_allele+alt_allele) != 0.5 & sum/length != 1))/nrow(subset(data3_map3,ref_allele/(ref_allele+alt_allele) != 0.5));

perfect <- subset(data3_map3,sum/length == 1);
imperfect <- subset(data3_map3,sum/length != 1);
perfect_hist <- hist(perfect$ref_allele/(perfect$ref_allele+perfect$alt_allele),breaks=seq(0,1,0.05));
imperfect_hist <- hist(imperfect$ref_allele/(imperfect$ref_allele+imperfect$alt_allele),breaks=seq(0,1,0.05));
mat <- cbind(perfect_hist$counts/nrow(perfect),imperfect_hist$counts/nrow(imperfect));

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_3mm_barplot_by_mapp.pdf");
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),ylim=c(0,1),xlab="",ylab="",main="",col=c("white","grey"));
dev.off();


data0_map0 <- merge(data0,map0,by.x=c("chr","pos"),by.y=c("locus","position"));
write.table(data0_map0,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_SNPs_36b_0mm.txt",quote=F,sep="\t",row.names=F,col.names=F);
perfect <- subset(data0_map0,((sum.x/length.x)+(sum.y/length.y)) == 2);
write.table(perfect,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_SNPs_36b_0mm_perfect.txt",quote=F,sep="\t",row.names=F,col.names=F);
nrow(perfect);

nrow(subset(data0_map0,ref_allele.x/(ref_allele.x+ref_allele.y) != 0.5 & ((sum.x/length.x)+(sum.y/length.y)) 1= 2))/nrow(subset(data0_map0,ref_allele.x/(ref_allele.x+ref_allele.y) != 0.5));

perfect <- subset(data0_map0,((sum.x/length.x)+(sum.y/length.y)) == 2);
imperfect <- subset(data0_map0,((sum.x/length.x)+(sum.y/length.y)) != 2);
perfect_hist <- hist(perfect$ref_allele.x/(perfect$ref_allele.x+perfect$ref_allele.y),breaks=seq(0,1,0.05));
imperfect_hist <- hist(imperfect$ref_allele.x/(imperfect$ref_allele.x+imperfect$ref_allele.y),breaks=seq(0,1,0.05));
mat <- cbind(perfect_hist$counts/nrow(perfect),imperfect_hist$counts/nrow(imperfect));

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_0mm_barplot_by_mapp.pdf");
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),ylim=c(0,1),xlab="",ylab="",main="",col=c("white","grey"));
dev.off();

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



######### integrate mappability and indels ############

library(VennDiagram);

map1 <- read.table("/Users/kraigrs/Wittkopp/Graze/berlin-updated-exonic-regions.SNPs.l36_m1.mappability.txt",header=TRUE,sep="\t");
map2 <- read.table("/Users/kraigrs/Wittkopp/Graze/berlin-updated-exonic-regions.SNPs.l36_m2.mappability.txt",header=TRUE,sep="\t");
map3 <- read.table("/Users/kraigrs/Wittkopp/Graze/berlin-updated-exonic-regions.SNPs.l36_m3.mappability.txt",header=TRUE,sep="\t");

berlin_map <- read.table("/Users/kraigrs/Wittkopp/Graze/berlin-updated-exonic-regions_fsa_masked.SNPs.mappability.txt",header=TRUE,sep="\t");
c1674_map <- read.table("/Users/kraigrs/Wittkopp/Graze/c1674-updated-exonic-regions_fsa_masked.SNPs.mappability.txt",header=TRUE,sep="\t");
map0 <- merge(berlin_map,c1674_map,by.x=c("locus","position"),by.y=c("locus","position"));

temp <- read.table("/Users/kraigrs/Wittkopp/Graze/berlin_c1674_fsa.indels.txt",header=TRUE,sep="\t");
key <- read.table("/Users/kraigrs/Wittkopp/Graze/key",header=TRUE,sep="\t");

indels <- merge(temp,key,by.x=c("chr","pos"),by.y=c("locus","fsa"));

# 1 mismatch
data1_map1 <- merge(data1,map1,by.x=c("chr","pos"),by.y=c("locus","position"));

data1_map1_indels <- merge(data1_map1,indels,by.x=c("chr","pos"),by.y=c("chr","mel"));

#nrow(subset(data1_map1_indels,neighbors < 1))/nrow(data1_map1_indels);

nrow(subset(data1_map1_indels,neighbors >= 1));
nrow(subset(data1_map1_indels,sum/length != 1));
nrow(subset(data1_map1_indels,ref_indel + alt_indel != 0));
nrow(subset(data1_map1_indels,neighbors >= 1 & sum/length != 1));
nrow(subset(data1_map1_indels,neighbors >= 1 & ref_indel + alt_indel != 0));
nrow(subset(data1_map1_indels,sum/length != 1 & ref_indel + alt_indel != 0));
nrow(subset(data1_map1_indels,neighbors >= 1 & sum/length != 1 & ref_indel + alt_indel != 0));

area1 <- nrow(subset(data1_map1_indels,neighbors >= 1));
area2 <- nrow(subset(data1_map1_indels,sum/length != 1));
area3 <- nrow(subset(data1_map1_indels,ref_indel + alt_indel != 0));
n12 <- nrow(subset(data1_map1_indels,neighbors >= 1 & sum/length != 1));
n13 <- nrow(subset(data1_map1_indels,neighbors >= 1 & ref_indel + alt_indel != 0));
n23 <- nrow(subset(data1_map1_indels,sum/length != 1 & ref_indel + alt_indel != 0));
n123 <- nrow(subset(data1_map1_indels,neighbors >= 1 & sum/length != 1 & ref_indel + alt_indel != 0));

venn.plot <- draw.triple.venn(area1,area2,area3,n12,n23,n13,n123,
	category = c("Biased", "Imperfect mappability", "Indels"),
	scaled = TRUE,
	fill = c("cyan", "magenta", "yellow"),
	lty = "blank",
	cex = 2,
	cat.cex = 2,
	cat.col = "black"
);

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/melsim_1mm_venn.pdf");
grid.draw(venn.plot);
dev.off();


nrow(subset(data1_map1_indels,neighbors < 1 & sum/length == 1 & ref_indel + alt_indel == 0));
nrow(subset(data1_map1_indels,neighbors < 1 & sum/length == 1 & ref_indel + alt_indel == 0 & ref_allele == alt_allele));
nrow(subset(data1_map1_indels,neighbors < 1 & sum/length == 1 & ref_indel + alt_indel == 0 & ref_allele != alt_allele));

nrow(data1_map1);
nrow(subset(data1_map1,ref_allele == alt_allele));

nrow(subset(data1_map1,neighbors < 1));
nrow(subset(data1_map1,neighbors < 1 & ref_allele == alt_allele));


nrow(subset(data1_map1_indels,neighbors < 1 & sum/length == 1));
nrow(subset(data1_map1_indels,neighbors < 1 & sum/length == 1 & ref_allele == alt_allele));
nrow(subset(data1_map1_indels,neighbors < 1 & sum/length == 1 & ref_allele != alt_allele));

nrow(subset(data1_map1_indels,neighbors < 1 & sum/length == 1 & ref_indel + alt_indel == 0));
nrow(subset(data1_map1_indels,neighbors < 1 & sum/length == 1 & ref_indel + alt_indel == 0 & ref_allele == alt_allele));
nrow(subset(data1_map1_indels,neighbors < 1 & sum/length == 1 & ref_indel + alt_indel == 0 & ref_allele != alt_allele));

write.table(subset(data1_map1_indels,neighbors < 1 & sum/length == 1 & ref_indel + alt_indel == 0),file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_SNPs_36b_1mm_no_indel.txt",quote=F,sep="\t",row.names=F,col.names=F);


AI <- subset(data1_map1_indels,ref_allele != alt_allele);
nrow(subset(AI,neighbors < 1 & sum/length != 1))/nrow(subset(AI,neighbors < 1));
nrow(subset(AI,neighbors < 1 & sum/length == 1 & ref_indel + alt_indel >= 1))/nrow(subset(AI,neighbors < 1 & sum/length == 1));

nrow(subset(AI,neighbors >= 1))/nrow(AI);
nrow(subset(AI,sum/length != 1))/nrow(AI);
nrow(subset(AI,ref_indel + alt_indel >= 1))/nrow(AI);


unbiased <- subset(data1_map1_indels,neighbors<1);
nrow(unbiased);
nrow(subset(unbiased,ref_allele != alt_allele));
nrow(subset(unbiased,ref_allele != alt_allele & sum/length != 1));
nrow(subset(unbiased,ref_allele != alt_allele & ref_indel + alt_indel != 0));


perfect <- subset(unbiased,sum/length == 1);
imperfect <- subset(unbiased,sum/length != 1);
perfect_hist <- hist(perfect$ref_allele/(perfect$ref_allele+perfect$alt_allele),breaks=seq(0,1,0.05));
imperfect_hist <- hist(imperfect$ref_allele/(imperfect$ref_allele+imperfect$alt_allele),breaks=seq(0,1,0.05));
mat <- cbind(perfect_hist$counts/nrow(perfect),imperfect_hist$counts/nrow(imperfect));

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_1mm_barplot_by_mapp.pdf");
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),ylim=c(0,1),xlab="",ylab="",main="",col=c("white","grey"));
dev.off();

no_indel <- subset(perfect,ref_indel + alt_indel == 0);
indel <- subset(perfect,ref_indel + alt_indel != 0);
no_indel_hist <- hist(no_indel$ref_allele/(no_indel$ref_allele+no_indel$alt_allele),breaks=seq(0,1,0.05));
indel_hist <- hist(indel$ref_allele/(indel$ref_allele+indel$alt_allele),breaks=seq(0,1,0.05));
mat <- cbind(no_indel_hist$counts/nrow(no_indel),indel_hist$counts/nrow(indel));

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_1mm_barplot_by_indel.pdf");
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),ylim=c(0,1),xlab="",ylab="",main="",col=c("white","grey"));
dev.off();

# 2 mismatches
data2_map2 <- merge(data2,map2,by.x=c("chr","pos"),by.y=c("locus","position"));
data2_map2_indels <- merge(data2_map2,indels,by.x=c("chr","pos"),by.y=c("chr","mel"));

nrow(subset(data2_map2_indels,neighbors >= 2));
nrow(subset(data2_map2_indels,sum/length != 1));
nrow(subset(data2_map2_indels,ref_indel + alt_indel != 0));
nrow(subset(data2_map2_indels,neighbors >= 2 & sum/length != 1));
nrow(subset(data2_map2_indels,neighbors >= 2 & ref_indel + alt_indel != 0));
nrow(subset(data2_map2_indels,sum/length != 1 & ref_indel + alt_indel != 0));
nrow(subset(data2_map2_indels,neighbors >= 2 & sum/length != 1 & ref_indel + alt_indel != 0));

nrow(data2_map2_indels);
nrow(subset(data2_map2_indels,neighbors < 2));
nrow(subset(data2_map2_indels,sum/length == 1));
nrow(subset(data2_map2_indels,ref_indel + alt_indel == 0));
nrow(subset(data2_map2_indels,neighbors < 2 & sum/length == 1));
nrow(subset(data2_map2_indels,neighbors < 2 & ref_indel + alt_indel == 0));
nrow(subset(data2_map2_indels,sum/length == 1 & ref_indel + alt_indel == 0));
nrow(subset(data2_map2_indels,neighbors < 2 & sum/length == 1 & ref_indel + alt_indel == 0));

area1 <- nrow(subset(data2_map2_indels,neighbors >= 2));
area2 <- nrow(subset(data2_map2_indels,sum/length != 1));
area3 <- nrow(subset(data2_map2_indels,ref_indel + alt_indel != 0));
n12 <- nrow(subset(data2_map2_indels,neighbors >= 2 & sum/length != 1));
n13 <- nrow(subset(data2_map2_indels,neighbors >= 2 & ref_indel + alt_indel != 0));
n23 <- nrow(subset(data2_map2_indels,sum/length != 1 & ref_indel + alt_indel != 0));
n123 <- nrow(subset(data2_map2_indels,neighbors >= 2 & sum/length != 1 & ref_indel + alt_indel != 0));

venn.plot <- draw.triple.venn(area1,area2,area3,n12,n23,n13,n123,
	scaled = TRUE,
	fill = c("cyan", "magenta", "yellow"),
	lty = "blank",
	cex = 2,
	cat.cex = 2,
	cat.col = "black"
);

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/melsim_2mm_venn.pdf");
grid.draw(venn.plot);
dev.off();

nrow(data2_map2);
nrow(subset(data2_map2,ref_allele == alt_allele));

nrow(subset(data2_map2,neighbors < 2));
nrow(subset(data2_map2,neighbors < 2 & ref_allele == alt_allele));

nrow(subset(data2_map2,ref_allele == alt_allele & sum/length == 1 & neighbors < 2))/nrow(subset(data2_map2,sum/length == 1 & neighbors < 2));

nrow(subset(data2_map2_indels,neighbors < 2));

nrow(subset(data2_map2_indels,neighbors < 2 & sum/length == 1));
nrow(subset(data2_map2_indels,neighbors < 2 & sum/length == 1 & ref_allele == alt_allele));
nrow(subset(data2_map2_indels,neighbors < 2 & sum/length == 1 & ref_allele != alt_allele));

nrow(subset(data2_map2_indels,neighbors < 2 & sum/length == 1 & ref_indel + alt_indel == 0));
nrow(subset(data2_map2_indels,neighbors < 2 & sum/length == 1 & ref_indel + alt_indel == 0 & ref_allele == alt_allele));
nrow(subset(data2_map2_indels,neighbors < 2 & sum/length == 1 & ref_indel + alt_indel == 0 & ref_allele != alt_allele));


write.table(subset(data2_map2_indels,neighbors < 2 & sum/length == 1 & ref_indel + alt_indel == 0),file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_SNPs_36b_2mm_no_indel.txt",quote=F,sep="\t",row.names=F,col.names=F);

AI <- subset(data2_map2_indels,ref_allele != alt_allele);
nrow(subset(AI,neighbors < 2 & sum/length != 1))/nrow(subset(AI,neighbors < 2));
nrow(subset(AI,neighbors < 2 & sum/length == 1 & ref_indel + alt_indel >= 1))/nrow(subset(AI,neighbors < 2 & sum/length == 1));

unbiased <- subset(data2_map2_indels,neighbors<2);
nrow(unbiased);
nrow(subset(unbiased,ref_allele != alt_allele));
nrow(subset(unbiased,ref_allele != alt_allele & sum/length != 1));
nrow(subset(unbiased,ref_allele != alt_allele & ref_indel + alt_indel != 0));


perfect <- subset(unbiased,sum/length == 1);
imperfect <- subset(unbiased,sum/length != 1);
perfect_hist <- hist(perfect$ref_allele/(perfect$ref_allele+perfect$alt_allele),breaks=seq(0,1,0.05));
imperfect_hist <- hist(imperfect$ref_allele/(imperfect$ref_allele+imperfect$alt_allele),breaks=seq(0,1,0.05));
mat <- cbind(perfect_hist$counts/nrow(perfect),imperfect_hist$counts/nrow(imperfect));

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_2mm_barplot_by_mapp.pdf");
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),ylim=c(0,1),xlab="",ylab="",main="",col=c("white","grey"));
dev.off();

no_indel <- subset(perfect,ref_indel + alt_indel == 0);
indel <- subset(perfect,ref_indel + alt_indel != 0);
no_indel_hist <- hist(no_indel$ref_allele/(no_indel$ref_allele+no_indel$alt_allele),breaks=seq(0,1,0.05));
indel_hist <- hist(indel$ref_allele/(indel$ref_allele+indel$alt_allele),breaks=seq(0,1,0.05));
mat <- cbind(no_indel_hist$counts/nrow(no_indel),indel_hist$counts/nrow(indel));

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_2mm_barplot_by_indel.pdf");
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),ylim=c(0,1),xlab="",ylab="",main="",col=c("white","grey"));
dev.off();



# 3 mismatches
data3_map3 <- merge(data3,map3,by.x=c("chr","pos"),by.y=c("locus","position"));
data3_map3_indels <- merge(data3_map3,indels,by.x=c("chr","pos"),by.y=c("chr","mel"));

nrow(subset(data3_map3_indels,neighbors >= 3));
nrow(subset(data3_map3_indels,sum/length != 1));
nrow(subset(data3_map3_indels,ref_indel + alt_indel != 0));
nrow(subset(data3_map3_indels,neighbors >= 3 & sum/length != 1));
nrow(subset(data3_map3_indels,neighbors >= 3 & ref_indel + alt_indel != 0));
nrow(subset(data3_map3_indels,sum/length != 1 & ref_indel + alt_indel != 0));
nrow(subset(data3_map3_indels,neighbors >= 3 & sum/length != 1 & ref_indel + alt_indel != 0));


nrow(data3_map3_indels);
nrow(subset(data3_map3_indels,neighbors < 3));
nrow(subset(data3_map3_indels,sum/length == 1));
nrow(subset(data3_map3_indels,ref_indel + alt_indel == 0));
nrow(subset(data3_map3_indels,neighbors < 3 & sum/length == 1));
nrow(subset(data3_map3_indels,neighbors < 3 & ref_indel + alt_indel == 0));
nrow(subset(data3_map3_indels,sum/length == 1 & ref_indel + alt_indel == 0));
nrow(subset(data3_map3_indels,neighbors < 3 & sum/length == 1 & ref_indel + alt_indel == 0));

area1 <- nrow(subset(data3_map3_indels,neighbors >= 3));
area2 <- nrow(subset(data3_map3_indels,sum/length != 1));
area3 <- nrow(subset(data3_map3_indels,ref_indel + alt_indel != 0));
n12 <- nrow(subset(data3_map3_indels,neighbors >= 3 & sum/length != 1));
n13 <- nrow(subset(data3_map3_indels,neighbors >= 3 & ref_indel + alt_indel != 0));
n23 <- nrow(subset(data3_map3_indels,sum/length != 1 & ref_indel + alt_indel != 0));
n123 <- nrow(subset(data3_map3_indels,neighbors >= 3 & sum/length != 1 & ref_indel + alt_indel != 0));

venn.plot <- draw.triple.venn(area1,area2,area3,n12,n23,n13,n123,
	scaled = TRUE,
	fill = c("cyan", "magenta", "yellow"),
	lty = "blank",
	cex = 2,
	cat.cex = 2,
	cat.col = "black"
);

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/melsim_3mm_venn.pdf");
grid.draw(venn.plot);
dev.off();


nrow(data3_map3);
nrow(subset(data3_map3,ref_allele == alt_allele));

nrow(subset(data3_map3,neighbors < 3));
nrow(subset(data3_map3,neighbors < 3 & ref_allele == alt_allele));

nrow(subset(data3_map3,ref_allele == alt_allele & sum/length == 1 & neighbors < 3))/nrow(subset(data3_map3,sum/length == 1 & neighbors < 3));

nrow(subset(data3_map3_indels,neighbors < 3));

nrow(subset(data3_map3_indels,neighbors < 3 & sum/length == 1));
nrow(subset(data3_map3_indels,neighbors < 3 & sum/length == 1 & ref_allele == alt_allele));
nrow(subset(data3_map3_indels,neighbors < 3 & sum/length == 1 & ref_allele != alt_allele));

nrow(subset(data3_map3_indels,neighbors < 3 & sum/length == 1 & ref_indel + alt_indel == 0));
nrow(subset(data3_map3_indels,neighbors < 3 & sum/length == 1 & ref_indel + alt_indel == 0 & ref_allele == alt_allele));
nrow(subset(data3_map3_indels,neighbors < 3 & sum/length == 1 & ref_indel + alt_indel == 0 & ref_allele != alt_allele));

write.table(subset(data3_map3_indels,neighbors < 3 & sum/length == 1 & ref_indel + alt_indel == 0),file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_SNPs_36b_3mm_no_indel.txt",quote=F,sep="\t",row.names=F,col.names=F);


AI <- subset(data3_map3_indels,ref_allele != alt_allele);
nrow(subset(AI,neighbors < 3 & sum/length != 1))/nrow(subset(AI,neighbors < 3));
nrow(subset(AI,neighbors < 3 & sum/length == 1 & ref_indel + alt_indel >= 1))/nrow(subset(AI,neighbors < 3 & sum/length == 1));

unbiased <- subset(data3_map3_indels,neighbors<3);
nrow(unbiased);
nrow(subset(unbiased,ref_allele != alt_allele));
nrow(subset(unbiased,ref_allele != alt_allele & sum/length != 1));
nrow(subset(unbiased,ref_allele != alt_allele & ref_indel + alt_indel != 0));

perfect <- subset(unbiased,sum/length == 1);
imperfect <- subset(unbiased,sum/length != 1);
perfect_hist <- hist(perfect$ref_allele/(perfect$ref_allele+perfect$alt_allele),breaks=seq(0,1,0.05));
imperfect_hist <- hist(imperfect$ref_allele/(imperfect$ref_allele+imperfect$alt_allele),breaks=seq(0,1,0.05));
mat <- cbind(perfect_hist$counts/nrow(perfect),imperfect_hist$counts/nrow(imperfect));

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_3mm_barplot_by_mapp.pdf");
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),ylim=c(0,1),xlab="",ylab="",main="",col=c("white","grey"));
dev.off();

no_indel <- subset(perfect,ref_indel + alt_indel == 0);
indel <- subset(perfect,ref_indel + alt_indel != 0);
no_indel_hist <- hist(no_indel$ref_allele/(no_indel$ref_allele+no_indel$alt_allele),breaks=seq(0,1,0.05));
indel_hist <- hist(indel$ref_allele/(indel$ref_allele+indel$alt_allele),breaks=seq(0,1,0.05));
mat <- cbind(no_indel_hist$counts/nrow(no_indel),indel_hist$counts/nrow(indel));

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_3mm_barplot_by_indel.pdf");
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),ylim=c(0,1),xlab="",ylab="",main="",col=c("white","grey"));
dev.off();


# 0 mismatches
data0_map0 <- merge(data0,map0,by.x=c("chr","pos"),by.y=c("locus","position"));
data0_map0_indels <- merge(data0_map0,indels,by.x=c("chr","pos"),by.y=c("chr","pos"));

nrow(data0_map0_indels);
nrow(subset(data0_map0_indels,((sum.x/length.x)+(sum.y/length.y)) != 2));


nrow(data0_map0);
nrow(subset(data0_map0,ref_allele.x == ref_allele.y));


nrow(subset(data0_map0,ref_allele.x == ref_allele.y & ((sum.x/length.x)+(sum.y/length.y)) == 2))/nrow(subset(data0_map0,((sum.x/length.x)+(sum.y/length.y)) == 2 ));

nrow(data0_map0_indels);

nrow(subset(data0_map0_indels,((sum.x/length.x)+(sum.y/length.y)) == 2));
nrow(subset(data0_map0_indels,((sum.x/length.x)+(sum.y/length.y)) == 2 & ref_allele.x == ref_allele.y));
nrow(subset(data0_map0_indels,((sum.x/length.x)+(sum.y/length.y)) == 2 & ref_allele.x != ref_allele.y));

nrow(subset(data0_map0_indels,((sum.x/length.x)+(sum.y/length.y)) == 2 & ref_indel + alt_indel == 0));
nrow(subset(data0_map0_indels,((sum.x/length.x)+(sum.y/length.y)) == 2 & ref_indel + alt_indel == 0 & ref_allele.x == ref_allele.y));
nrow(subset(data0_map0_indels,((sum.x/length.x)+(sum.y/length.y)) == 2 & ref_indel + alt_indel == 0 & ref_allele.x != ref_allele.y));


write.table(subset(data0_map0_indels,((sum.x/length.x)+(sum.y/length.y)) == 2 & ref_indel + alt_indel == 0),file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_SNPs_36b_0mm_no_indel.txt",quote=F,sep="\t",row.names=F,col.names=F);


AI <- subset(data0_map0_indels,ref_allele.x != ref_allele.y);
nrow(subset(AI,((sum.x/length.x)+(sum.y/length.y)) != 2))/nrow(AI);
nrow(subset(AI,((sum.x/length.x)+(sum.y/length.y)) == 2 & ref_indel + alt_indel >= 1))/nrow(subset(AI,((sum.x/length.x)+(sum.y/length.y)) == 2));

unbiased <- data0_map0_indels;
nrow(unbiased);
nrow(subset(unbiased,ref_allele.x != ref_allele.y));
nrow(subset(unbiased,ref_allele.x != ref_allele.y & ((sum.x/length.x)+(sum.y/length.y)) != 2));
nrow(subset(unbiased,ref_allele.x != ref_allele.y & ref_indel + alt_indel != 0));

perfect <- subset(unbiased,((sum.x/length.x)+(sum.y/length.y)) == 2);
imperfect <- subset(unbiased,((sum.x/length.x)+(sum.y/length.y)) != 2);
perfect_hist <- hist(perfect$ref_allele.x/(perfect$ref_allele.x+perfect$ref_allele.y),breaks=seq(0,1,0.05));
imperfect_hist <- hist(imperfect$ref_allele.x/(imperfect$ref_allele.x+imperfect$ref_allele.y),breaks=seq(0,1,0.05));
mat <- cbind(perfect_hist$counts/nrow(perfect),imperfect_hist$counts/nrow(imperfect));

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_0mm_barplot_by_mapp.pdf");
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),ylim=c(0,1),xlab="",ylab="",main="",col=c("white","grey"));
dev.off();

no_indel <- subset(perfect,ref_indel + alt_indel == 0);
indel <- subset(perfect,ref_indel + alt_indel != 0);
no_indel_hist <- hist(no_indel$ref_allele.x/(no_indel$ref_allele.x+no_indel$ref_allele.y),breaks=seq(0,1,0.05));
indel_hist <- hist(indel$ref_allele.x/(indel$ref_allele.x+indel$ref_allele.y),breaks=seq(0,1,0.05));
mat <- cbind(no_indel_hist$counts/nrow(no_indel),indel_hist$counts/nrow(indel));

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_0mm_barplot_by_indel.pdf");
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),ylim=c(0,1),xlab="",ylab="",main="",col=c("white","grey"));
dev.off();

######### make plots comparing approaches ############

key <- read.table("/Users/kraigrs/Wittkopp/Graze/key",header=TRUE,sep="\t");

perfect0 <- subset(data0_tests_map0,((sum.x/length.x)+(sum.y/length.y)) == 2);

temp1 <- subset(data1_tests_map1,neighbors<1);
perfect1 <- subset(temp1,sum/length == 1);

temp2 <- subset(data2_tests_map2,neighbors<2);
perfect2 <- subset(temp2,sum/length == 1);

#temp3 <- subset(data3_tests_map3,neighbors<3);
temp3 <- data3_tests_map3;
perfect3 <- subset(temp3,sum/length == 1);


x <- merge(perfect0,key,by.x=c("chr","pos"),by.y=c("locus","fsa"));
y <- merge(perfect3,key,by.x=c("chr","pos"),by.y=c("locus","mel"));
compare <- merge(x,y,by.x=c("chr","mel"),by.y=c("chr","pos"));
cor(compare$berlin_ref_allele/(compare$berlin_ref_allele+compare$c1674_alt_allele) , compare$ref_allele/(compare$ref_allele+compare$alt_allele))^2;

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

pdf(file="/Users/kraigrs/Desktop/new_plots2/compare_real_data_3mm_biased.pdf");

plot( compare$berlin_ref_allele/(compare$berlin_ref_allele+compare$c1674_alt_allele) , compare$ref_allele/(compare$ref_allele+compare$alt_allele) ,xlim=c(0,1),ylim=c(0,1),pch=19,cex=0.5,col=rgb(0,0,0,0.05),xlab="",ylab="",main="");

dev.off();


pdf(file="/Users/kraigrs/Desktop/plots/individual/pie_real_data_3mm_biased.pdf");

n_n <- subset(compare,qvals.x >= 0.05 & qvals.y >= 0.05);
s_ns <- subset(compare,qvals.x < 0.05 & qvals.y >= 0.05);
ns_s <- subset(compare,qvals.x >= 0.05 & qvals.y < 0.05);
s_s <- subset(compare,qvals.x < 0.05 & qvals.y < 0.05);

pie(c(nrow(n_n),nrow(s_ns),nrow(ns_s),nrow(s_s)),col=c("grey","red","blue","purple"),labels="");

dev.off();


plot( compare$berlin_ref_allele/(compare$berlin_ref_allele+compare$c1674_alt_allele) , compare$ref_allele/(compare$ref_allele+compare$alt_allele) ,xlim=c(0,1),ylim=c(0,1),type="n",pch=19,cex=0.5,col=rgb(0,0,0,0.05),xlab="",ylab="",main="");

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




