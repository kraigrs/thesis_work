data1 <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v1_m1.SNPs.txt",header=TRUE,sep="\t");

data2 <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v2_m1.SNPs.txt",header=TRUE,sep="\t");

data3 <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v3_m1.SNPs.txt",header=TRUE,sep="\t");

data0 <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.bowtie_v0_m1.SNPs.txt",header=TRUE,sep="\t");
DGRP <- read.table("/Users/kraigrs/Wittkopp/DGRP/DGRP_line_40_SNPs_const.txt",header=FALSE,sep="\t");
data0 <- merge(data0,DGRP,by.x=c("chr","pos"),by.y=c("V1","V3"));

data1_alt <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_alt_line_40.bowtie_v1_m1.SNPs.txt",header=TRUE,sep="\t");

relevant <- read.table("/Users/kraigrs/Wittkopp/Simulations/sim_regs.50b_windows_SNPs.neighbors.txt",header=TRUE,sep="\t");

summary(relevant$relevant);

data1 <- merge(data1,relevant,by.x=c("chr","pos"),by.y=c("chr","pos"));
data2 <- merge(data2,relevant,by.x=c("chr","pos"),by.y=c("chr","pos"));
data3 <- merge(data3,relevant,by.x=c("chr","pos"),by.y=c("chr","pos"));
data0 <- merge(data0,relevant,by.x=c("chr","pos"),by.y=c("chr","pos"));
data1_alt <- merge(data1_alt,relevant,by.x=c("chr","pos"),by.y=c("chr","pos"));

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_v1_box.pdf");
boxplot(ref_allele/(ref_allele+alt_allele) ~ relevant,data = data1, varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,11),ylim=c(0,1),pars = list(outpch=19,cex=0.2));
abline(h=0.5,lty=2,col="red");
dev.off();

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_v2_box.pdf");
boxplot(ref_allele/(ref_allele+alt_allele) ~ relevant,data = data2, varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,11),ylim=c(0,1),pars = list(outpch=19,cex=0.2));
abline(h=0.5,lty=2,col="red");
dev.off();

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_v3_box.pdf");
boxplot(ref_allele/(ref_allele+alt_allele) ~ relevant,data = data3, varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,11),ylim=c(0,1),pars = list(outpch=19,cex=0.2));
abline(h=0.5,lty=2,col="red");
dev.off();


pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_v0_box.pdf");
boxplot(dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele) ~ relevant,data = data0, varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,11),ylim=c(0,1),pars = list(outpch=19,cex=0.2));
abline(h=0.5,lty=2,col="red");
dev.off();

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_alt_v1_box.pdf");
boxplot(alt_allele/(ref_allele+alt_allele) ~ relevant,data = data1_alt, varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,11),ylim=c(0,1),pars = list(outpch=19,cex=0.2));
abline(h=0.5,lty=2,col="red");
dev.off();

### pie charts ###
pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_v1_pie.pdf");
pie(c(nrow(subset(data1,ref_allele/(ref_allele+alt_allele)!=0.5)),nrow(subset(data1,ref_allele/(ref_allele+alt_allele)==0.5))),col=c("gray","white"),labels="");
dev.off();

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_v2_pie.pdf");
pie(c(nrow(subset(data2,ref_allele/(ref_allele+alt_allele)!=0.5)),nrow(subset(data2,ref_allele/(ref_allele+alt_allele)==0.5))),col=c("gray","white"),labels="");
dev.off();

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_v3_pie.pdf");
pie(c(nrow(subset(data3,ref_allele/(ref_allele+alt_allele)!=0.5)),nrow(subset(data3,ref_allele/(ref_allele+alt_allele)==0.5))),col=c("gray","white"),labels="");
dev.off();

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_v0_pie.pdf");
pie(c(nrow(subset(data0,dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele)!=0.5)),nrow(subset(data0,dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele)==0.5))),col=c("gray","white"),labels="");
dev.off();

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_alt_v1_pie.pdf");
pie(c(nrow(subset(data1_alt,ref_allele/(ref_allele+alt_allele)!=0.5)),nrow(subset(data1_alt,ref_allele/(ref_allele+alt_allele)==0.5))),col=c("gray","white"),labels="");
dev.off();

nrow(subset(data1,ref_allele != alt_allele))/nrow(data1);
nrow(subset(data1,ref_allele > alt_allele))/nrow(subset(data1,ref_allele != alt_allele));

nrow(subset(data1,ref_allele > alt_allele & relevant > 1))/nrow(subset(data1,relevant > 1));

nrow(subset(data1,ref_allele == alt_allele))/nrow(data1);
nrow(subset(data2,ref_allele == alt_allele))/nrow(data2);
nrow(subset(data3,ref_allele == alt_allele))/nrow(data3);
nrow(subset(data0,dm3_ref_ref_allele == dm3_alt_alt_allele))/nrow(data0);


###############
# mappability #
###############

library(VennDiagram);

# 1 mismatch

mappability_1mm <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/SNPs_line_40.dm3_ref.l50_m1.mappability.txt",header=TRUE,sep="\t");

data <- merge(data1,mappability_1mm,by.x=c("chr","pos"),by.y=c("chr","position"));

nrow(subset(data,relevant > 1));
nrow(subset(data,sum/length != 1));
nrow(subset(data,relevant > 1 & sum/length != 1));

nrow(subset(data,relevant <= 1));
nrow(subset(data,sum/length == 1));
nrow(subset(data,relevant <= 1 & sum/length == 1));

write.table(data,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_SNPs_1mm.txt",quote=F,sep="\t",row.names=F,col.names=F);

temp <- subset(data,relevant <= 1);
write.table(temp,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_SNPs_1mm_unbiased.txt",quote=F,sep="\t",row.names=F,col.names=F);

perfect <- subset(temp,sum/length == 1);
imperfect <- subset(temp,sum/length != 1);

write.table(perfect,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_SNPs_1mm_unbiased_perfect.txt",quote=F,sep="\t",row.names=F,col.names=F);

nrow(subset(temp,ref_allele/(ref_allele+alt_allele) != 0.5 & sum/length != 1))/nrow(subset(temp,ref_allele/(ref_allele+alt_allele) != 0.5))*100; # % unequal ASRA and imperfect mappability 
nrow(subset(temp,ref_allele/(ref_allele+alt_allele) == 0.5 & sum/length != 1))/nrow(subset(temp,ref_allele/(ref_allele+alt_allele) == 0.5))*100; # % equal ASRA and imperfect mappability

# proportion of differentiating sites with perfect mappability and equal allelic abundance
nrow(subset(temp,ref_allele/(ref_allele+alt_allele) == 0.5 & sum/length == 1))/nrow(subset(temp,sum/length == 1))*100;


perfect_hist <- hist(perfect$ref_allele/(perfect$ref_allele+perfect$alt_allele),breaks=seq(0,1,0.05));
imperfect_hist <- hist(imperfect$ref_allele/(imperfect$ref_allele+imperfect$alt_allele),breaks=seq(0,1,0.05));
mat <- cbind(perfect_hist$counts/nrow(perfect),imperfect_hist$counts/nrow(imperfect));

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_1mm_barplot_by_mapp.pdf");
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),ylim=c(0,1),xlab="",ylab="",main="",col=c("white","grey"));
dev.off();


# 2 mismatches

mappability_2mm <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/SNPs_line_40.dm3_ref.l50_m2.mappability.txt",header=TRUE,sep="\t");

data <- merge(data2,mappability_2mm,by.x=c("chr","pos"),by.y=c("chr","position"));

nrow(subset(data,relevant > 2));
nrow(subset(data,sum/length != 1));
nrow(subset(data,relevant > 2 & sum/length != 1));

nrow(subset(data,relevant <= 2));
nrow(subset(data,sum/length == 1));
nrow(subset(data,relevant <= 2 & sum/length == 1));

write.table(data,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_SNPs_2mm.txt",quote=F,sep="\t",row.names=F,col.names=F);

temp <- subset(data,relevant <= 2);
write.table(temp,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_SNPs_2mm_unbiased.txt",quote=F,sep="\t",row.names=F,col.names=F);

perfect <- subset(temp,sum/length == 1);
imperfect <- subset(temp,sum/length != 1);

write.table(perfect,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_SNPs_2mm_unbiased_perfect.txt",quote=F,sep="\t",row.names=F,col.names=F);

nrow(subset(temp,ref_allele/(ref_allele+alt_allele) != 0.5 & sum/length != 1))/nrow(subset(temp,ref_allele/(ref_allele+alt_allele) != 0.5))*100; # % unequal ASRA and imperfect mappability 
nrow(subset(temp,ref_allele/(ref_allele+alt_allele) == 0.5 & sum/length != 1))/nrow(subset(temp,ref_allele/(ref_allele+alt_allele) == 0.5))*100; # % equal ASRA and imperfect mappability

# proportion of differentiating sites with perfect mappability and equal allelic abundance
nrow(subset(temp,ref_allele/(ref_allele+alt_allele) != 0.5 & sum/length != 1))/nrow(subset(temp,sum/length != 1))*100;


perfect_hist <- hist(perfect$ref_allele/(perfect$ref_allele+perfect$alt_allele),breaks=seq(0,1,0.05));
imperfect_hist <- hist(imperfect$ref_allele/(imperfect$ref_allele+imperfect$alt_allele),breaks=seq(0,1,0.05));
mat <- cbind(perfect_hist$counts/nrow(perfect),imperfect_hist$counts/nrow(imperfect));

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_2mm_barplot_by_mapp.pdf");
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),ylim=c(0,1),xlab="",ylab="",main="",col=c("white","grey"));
dev.off();


# 3 mismatches

mappability_3mm <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/SNPs_line_40.dm3_ref.l50_m3.mappability.txt",header=TRUE,sep="\t");

data <- merge(data3,mappability_3mm,by.x=c("chr","pos"),by.y=c("chr","position"));

nrow(subset(data,relevant > 3));
nrow(subset(data,sum/length != 1));
nrow(subset(data,relevant > 3 & sum/length != 1));

nrow(subset(data,relevant <= 3));
nrow(subset(data,sum/length == 1));
nrow(subset(data,relevant <= 3 & sum/length == 1));

write.table(data,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_SNPs_3mm.txt",quote=F,sep="\t",row.names=F,col.names=F);

temp <- subset(data,relevant <= 3);
write.table(temp,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_SNPs_3mm_unbiased.txt",quote=F,sep="\t",row.names=F,col.names=F);

perfect <- subset(temp,sum/length == 1);
imperfect <- subset(temp,sum/length != 1);

write.table(perfect,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_SNPs_3mm_unbiased_perfect.txt",quote=F,sep="\t",row.names=F,col.names=F);

nrow(subset(temp,ref_allele/(ref_allele+alt_allele) != 0.5 & sum/length != 1))/nrow(subset(temp,ref_allele/(ref_allele+alt_allele) != 0.5))*100; # % unequal ASRA and imperfect mappability 
nrow(subset(temp,ref_allele/(ref_allele+alt_allele) == 0.5 & sum/length != 1))/nrow(subset(temp,ref_allele/(ref_allele+alt_allele) == 0.5))*100; # % equal ASRA and imperfect mappability

# proportion of differentiating sites with perfect mappability and equal allelic abundance
nrow(subset(temp,ref_allele/(ref_allele+alt_allele) != 0.5 & sum/length != 1))/nrow(subset(temp,sum/length != 1))*100;


perfect_hist <- hist(perfect$ref_allele/(perfect$ref_allele+perfect$alt_allele),breaks=seq(0,1,0.05));
imperfect_hist <- hist(imperfect$ref_allele/(imperfect$ref_allele+imperfect$alt_allele),breaks=seq(0,1,0.05));
mat <- cbind(perfect_hist$counts/nrow(perfect),imperfect_hist$counts/nrow(imperfect));

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_3mm_barplot_by_mapp.pdf");
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),ylim=c(0,1),xlab="",ylab="",main="",col=c("white","grey"));
dev.off();


# 0 mismatches

ref_mappability <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/SNPs_line_40.dm3_ref.l50_m0.mappability.txt",header=TRUE,sep="\t");
alt_mappability <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/SNPs_line_40.dm3_alt_line_40.l50_m0.mappability.txt",header=TRUE,sep="\t");
mappability_0mm <- merge(ref_mappability,alt_mappability,by.x=c("chr","position"),by.y=c("chr","position"));

data <- merge(data0,mappability_0mm,by.x=c("chr","pos"),by.y=c("chr","position"));

nrow(data);
nrow(subset(data,((sum.x/length.x)+(sum.y/length.y)) == 2));

nrow(subset(data,dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele) != 0.5 & ((sum.x/length.x)+(sum.y/length.y)) != 2))/nrow(subset(data,((sum.x/length.x)+(sum.y/length.y)) != 2))*100;

write.table(data,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_SNPs_0mm.txt",quote=F,sep="\t",row.names=F,col.names=F);

perfect <- subset(data,((sum.x/length.x)+(sum.y/length.y)) == 2);
imperfect <- subset(data,((sum.x/length.x)+(sum.y/length.y)) != 2);

write.table(perfect,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_SNPs_0mm_unbiased_perfect.txt",quote=F,sep="\t",row.names=F,col.names=F);

# proportion of differentiating sites with perfect mappability and equal allelic abundance
nrow(subset(data,dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele) == 0.5 & ((sum.x/length.x)+(sum.y/length.y)) == 2))/nrow(subset(data,((sum.x/length.x)+(sum.y/length.y)) == 2))*100;


perfect_hist <- hist(perfect$dm3_ref_ref_allele/(perfect$dm3_ref_ref_allele+perfect$dm3_alt_alt_allele),breaks=seq(0,1,0.05));
imperfect_hist <- hist(imperfect$dm3_ref_ref_allele/(imperfect$dm3_ref_ref_allele+imperfect$dm3_alt_alt_allele),breaks=seq(0,1,0.05));
mat <- cbind(perfect_hist$counts/nrow(perfect),imperfect_hist$counts/nrow(imperfect));

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_50b_0mm_barplot_by_mapp.pdf");
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),ylim=c(0,1),xlab="",ylab="",main="",col=c("white","grey"));
dev.off();
