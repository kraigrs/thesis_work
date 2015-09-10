data1 <- read.table("/Users/kraigrs/Wittkopp/Graze/simulation/berlin_c1674.tiled_36bp.berlin-updated-exonic-regions.bowtie_v1_m1.SNPs.txt",header=TRUE,sep="\t");

data2 <- read.table("/Users/kraigrs/Wittkopp/Graze/simulation/berlin_c1674.tiled_36bp.berlin-updated-exonic-regions.bowtie_v2_m1.SNPs.txt",header=TRUE,sep="\t");

data3 <- read.table("/Users/kraigrs/Wittkopp/Graze/simulation/berlin_c1674.tiled_36bp.berlin-updated-exonic-regions.bowtie_v3_m1.SNPs.txt",header=TRUE,sep="\t");

data0_ref <- read.table("/Users/kraigrs/Wittkopp/Graze/simulation/berlin_c1674.tiled_36bp.berlin-updated-exonic-regions_fsa_masked.bowtie_v0_m1.SNPs.txt",header=TRUE,sep="\t");

data0_alt <- read.table("/Users/kraigrs/Wittkopp/Graze/simulation/berlin_c1674.tiled_36bp.c1674-updated-exonic-regions_fsa_masked.bowtie_v0_m1.SNPs.txt",header=TRUE,sep="\t");

data0 <- merge(data0_ref,data0_alt,by.x=c("chr","pos"),by.y=c("chr","pos"));

relevant <- read.table("/Users/kraigrs/Wittkopp/Graze/berlin-updated-exonic-regions_fsa.36b_windows_SNPs.neighbors.txt",header=TRUE,sep="\t");

key <- read.table("/Users/kraigrs/Wittkopp/Graze/key",header=TRUE,sep="\t");

summary(relevant$relevant);

relevant <- merge(relevant,key,by.x=c("chr","pos"),by.y=c("locus","mel"));

data1 <- merge(data1,relevant,by.x=c("chr","pos"),by.y=c("chr","pos"));
data2 <- merge(data2,relevant,by.x=c("chr","pos"),by.y=c("chr","pos"));
data3 <- merge(data3,relevant,by.x=c("chr","pos"),by.y=c("chr","pos"));
data0 <- merge(data0,relevant,by.x=c("chr","pos"),by.y=c("chr","fsa"));
data0 <- data0[,c(1:12,14,15)];

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_36b_v1_box.pdf");
boxplot( ref_allele/(ref_allele+alt_allele) ~ relevant, data = data1,varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,24),ylim=c(0,1),pars = list(outpch=19,outcex=0.2));
abline(h=0.5,lty=2,col="red");
dev.off();

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_36b_v2_box.pdf");
boxplot(ref_allele/(ref_allele+alt_allele) ~ relevant, data = data2,varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,24),ylim=c(0,1),pars = list(outpch=19,outcex=0.2));
abline(h=0.5,lty=2,col="red");
dev.off();

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_36b_v3_box.pdf");
boxplot(ref_allele/(ref_allele+alt_allele) ~ relevant, data = data3,varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,24),ylim=c(0,1),pars = list(outpch=19,outcex=0.2));
abline(h=0.5,lty=2,col="red");
dev.off();

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_36b_v0_box.pdf");
boxplot(ref_allele.x/(ref_allele.x+ref_allele.y) ~ relevant, data = data0,varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,24),ylim=c(0,1),pars = list(outpch=19,outcex=0.2));
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


######### integrate mappability and indels ############

map1 <- read.table("/Users/kraigrs/Wittkopp/Graze/berlin-updated-exonic-regions.SNPs.l36_m1.mappability.txt",header=TRUE,sep="\t");
map2 <- read.table("/Users/kraigrs/Wittkopp/Graze/berlin-updated-exonic-regions.SNPs.l36_m2.mappability.txt",header=TRUE,sep="\t");
map3 <- read.table("/Users/kraigrs/Wittkopp/Graze/berlin-updated-exonic-regions.SNPs.l36_m3.mappability.txt",header=TRUE,sep="\t");

berlin_map <- read.table("/Users/kraigrs/Wittkopp/Graze/berlin-updated-exonic-regions_fsa_masked.SNPs.mappability.txt",header=TRUE,sep="\t");
c1674_map <- read.table("/Users/kraigrs/Wittkopp/Graze/c1674-updated-exonic-regions_fsa_masked.SNPs.mappability.txt",header=TRUE,sep="\t");
map0 <- merge(berlin_map,c1674_map,by.x=c("locus","position"),by.y=c("locus","position"));

temp <- read.table("/Users/kraigrs/Wittkopp/Graze/berlin_c1674_fsa.indels.txt",header=TRUE,sep="\t");

indels <- merge(temp,key,by.x=c("chr","pos"),by.y=c("locus","fsa"));

# 1 mismatch
data1_map1 <- merge(data1,map1,by.x=c("chr","pos"),by.y=c("locus","position"));

unbiased <- subset(data1_map1,relevant <= 1);
perfect <- subset(unbiased,sum/length == 1);
nrow(subset(perfect,ref_allele == alt_allele))/nrow(perfect)*100;


data1_map1_indels <- merge(data1_map1,indels,by.x=c("chr","pos"),by.y=c("chr","mel"));

nrow(subset(data1_map1_indels,relevant > 1));
nrow(subset(data1_map1_indels,sum/length != 1));
nrow(subset(data1_map1_indels,ref_indel + alt_indel != 0));
nrow(subset(data1_map1_indels,relevant > 1 & sum/length != 1));
nrow(subset(data1_map1_indels,relevant > 1 & ref_indel + alt_indel != 0));
nrow(subset(data1_map1_indels,sum/length != 1 & ref_indel + alt_indel != 0));
nrow(subset(data1_map1_indels,relevant > 1 & sum/length != 1 & ref_indel + alt_indel != 0));

write.table(data1_map1_indels,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_SNPs_36b_1mm.txt",quote=F,sep="\t",row.names=F,col.names=F);

write.table(subset(data1_map1_indels,relevant <= 1),file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_SNPs_36b_1mm_unbiased.txt",quote=F,sep="\t",row.names=F,col.names=F);

write.table(subset(data1_map1_indels,relevant <= 1 & sum/length == 1),file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_SNPs_36b_1mm_perfect.txt",quote=F,sep="\t",row.names=F,col.names=F);

write.table(subset(data1_map1_indels,relevant <= 1 & sum/length == 1 & ref_indel + alt_indel == 0),file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_SNPs_36b_1mm_no_indel.txt",quote=F,sep="\t",row.names=F,col.names=F);

unbiased <- subset(data1_map1_indels,relevant <= 1);

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

unbiased <- subset(data2_map2,relevant <= 2);
perfect <- subset(unbiased,sum/length == 1);
nrow(subset(perfect,ref_allele == alt_allele))/nrow(perfect)*100;

data2_map2_indels <- merge(data2_map2,indels,by.x=c("chr","pos"),by.y=c("chr","mel"));

nrow(data2_map2_indels);
nrow(subset(data2_map2_indels,relevant > 2));
nrow(subset(data2_map2_indels,sum/length != 1));
nrow(subset(data2_map2_indels,ref_indel + alt_indel != 0));
nrow(subset(data2_map2_indels,relevant > 2 & sum/length != 1));
nrow(subset(data2_map2_indels,relevant > 2 & ref_indel + alt_indel != 0));
nrow(subset(data2_map2_indels,sum/length != 1 & ref_indel + alt_indel != 0));
nrow(subset(data2_map2_indels,relevant > 2 & sum/length != 1 & ref_indel + alt_indel != 0));


write.table(data2_map2_indels,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_SNPs_36b_2mm.txt",quote=F,sep="\t",row.names=F,col.names=F);

write.table(subset(data2_map2_indels,relevant <= 2),file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_SNPs_36b_2mm_unbiased.txt",quote=F,sep="\t",row.names=F,col.names=F);

write.table(subset(data2_map2_indels,relevant <= 2 & sum/length == 1),file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_SNPs_36b_2mm_perfect.txt",quote=F,sep="\t",row.names=F,col.names=F);

write.table(subset(data2_map2_indels,relevant <= 2 & sum/length == 1 & ref_indel + alt_indel == 0),file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_SNPs_36b_2mm_no_indel.txt",quote=F,sep="\t",row.names=F,col.names=F);

unbiased <- subset(data2_map2_indels,relevant <= 2);

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

unbiased <- subset(data3_map3,relevant <= 3);
perfect <- subset(unbiased,sum/length == 1);
nrow(subset(perfect,ref_allele == alt_allele))/nrow(perfect)*100;

data3_map3_indels <- merge(data3_map3,indels,by.x=c("chr","pos"),by.y=c("chr","mel"));

nrow(subset(data3_map3_indels,relevant > 3));
nrow(subset(data3_map3_indels,sum/length != 1));
nrow(subset(data3_map3_indels,ref_indel + alt_indel != 0));
nrow(subset(data3_map3_indels,relevant > 3 & sum/length != 1));
nrow(subset(data3_map3_indels,relevant > 3 & ref_indel + alt_indel != 0));
nrow(subset(data3_map3_indels,sum/length != 1 & ref_indel + alt_indel != 0));
nrow(subset(data3_map3_indels,relevant > 3 & sum/length != 1 & ref_indel + alt_indel != 0));

write.table(data3_map3_indels,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_SNPs_36b_3mm.txt",quote=F,sep="\t",row.names=F,col.names=F);

write.table(subset(data3_map3_indels,relevant <= 3),file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_SNPs_36b_3mm_unbiased.txt",quote=F,sep="\t",row.names=F,col.names=F);

write.table(subset(data3_map3_indels,relevant <= 3 & sum/length == 1),file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_SNPs_36b_3mm_perfect.txt",quote=F,sep="\t",row.names=F,col.names=F);

write.table(subset(data3_map3_indels,relevant <= 3 & sum/length == 1 & ref_indel + alt_indel == 0),file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_SNPs_36b_3mm_no_indel.txt",quote=F,sep="\t",row.names=F,col.names=F);

unbiased <- subset(data3_map3_indels,relevant <= 3);

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
#data0 <- data0[,c(1:12,14,15)];

data0_map0 <- merge(data0,map0,by.x=c("chr","pos"),by.y=c("locus","position"));

perfect <- subset(data0_map0,((sum.x/length.x)+(sum.y/length.y)) == 2);
nrow(subset(perfect,ref_allele.x == ref_allele.y))/nrow(perfect)*100;

data0_map0_indels <- merge(data0_map0,indels,by.x=c("chr","pos"),by.y=c("chr","pos"));

nrow(data0_map0_indels);
nrow(subset(data0_map0_indels,((sum.x/length.x)+(sum.y/length.y)) != 2));


write.table(data0_map0_indels,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_SNPs_36b_0mm.txt",quote=F,sep="\t",row.names=F,col.names=F);

write.table(subset(data0_map0_indels,((sum.x/length.x)+(sum.y/length.y)) == 2),file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/real_SNPs_36b_0mm_perfect.txt",quote=F,sep="\t",row.names=F,col.names=F);


unbiased <- data0_map0_indels;

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
