data1 <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp36_error0_tiled_line_40.dm3_ref.bowtie_v1_m1.SNPs.txt",header=TRUE,sep="\t");

data2 <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp36_error0_tiled_line_40.dm3_ref.bowtie_v2_m1.SNPs.txt",header=TRUE,sep="\t");

data3 <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp36_error0_tiled_line_40.dm3_ref.bowtie_v3_m1.SNPs.txt",header=TRUE,sep="\t");

data0_ref <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp36_error0_tiled_line_40.dm3_ref.bowtie_v0_m1.SNPs.txt",header=TRUE,sep="\t");

data0_alt <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp36_error0_tiled_line_40.dm3_alt_line_40.bowtie_v0_m1.SNPs.txt",header=TRUE,sep="\t");

data0 <- merge(data0_ref,data0_alt,by.x=c("chr","pos"),by.y=c("chr","pos"));

data1_alt <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp36_error0_tiled_line_40.dm3_alt_line_40.bowtie_v1_m1.SNPs.txt",header=TRUE,sep="\t");

nrow(subset(data1,ref_allele != alt_allele))/nrow(data1);
nrow(subset(data1,ref_allele != alt_allele & ref_allele > alt_allele))/nrow(subset(data1,ref_allele != alt_allele));
nrow(subset(data1,neighbors > 0))/nrow(data1);
nrow(subset(data1,neighbors > 0 & ref_allele > alt_allele))/nrow(subset(data1,neighbors > 0));

nrow(subset(data2,ref_allele == alt_allele))/nrow(data2);
nrow(subset(data3,ref_allele == alt_allele))/nrow(data3);

nrow(subset(data0,ref_allele.x == ref_allele.y))/nrow(data0);

### pie charts ###
pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_36b_v1_pie.pdf");
pie(c(nrow(subset(data1,ref_allele/(ref_allele+alt_allele)!=0.5)),nrow(subset(data1,ref_allele/(ref_allele+alt_allele)==0.5))),col=c("gray","white"),labels="");
dev.off();

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_36b_v2_pie.pdf");
pie(c(nrow(subset(data2,ref_allele/(ref_allele+alt_allele)!=0.5)),nrow(subset(data2,ref_allele/(ref_allele+alt_allele)==0.5))),col=c("gray","white"),labels="");
dev.off();

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_36b_v3_pie.pdf");
pie(c(nrow(subset(data3,ref_allele/(ref_allele+alt_allele)!=0.5)),nrow(subset(data3,ref_allele/(ref_allele+alt_allele)==0.5))),col=c("gray","white"),labels="");
dev.off();

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_36b_v0_pie.pdf");
pie(c(nrow(subset(data0,ref_allele.x/(ref_allele.x+ref_allele.y)!=0.5)),nrow(subset(data0,ref_allele.x/(ref_allele.x+ref_allele.y)==0.5))),col=c("gray","white"),labels="");
dev.off();

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_36b_alt_v1_pie.pdf");
pie(c(nrow(subset(data1_alt,ref_allele/(ref_allele+alt_allele)!=0.5)),nrow(subset(data1_alt,ref_allele/(ref_allele+alt_allele)==0.5))),col=c("gray","white"),labels="");
dev.off();

nrow(subset(data1,ref_allele/(ref_allele+alt_allele)!=0.5))/nrow(subset(data1,ref_allele/(ref_allele+alt_allele)==0.5));
nrow(subset(data2,ref_allele/(ref_allele+alt_allele)!=0.5))/nrow(subset(data1,ref_allele/(ref_allele+alt_allele)==0.5));
nrow(subset(data3,ref_allele/(ref_allele+alt_allele)!=0.5))/nrow(subset(data1,ref_allele/(ref_allele+alt_allele)==0.5));
nrow(subset(data0,ref_allele.x/(ref_allele.x+ref_allele.y)!=0.5))/nrow(subset(data0,ref_allele.x/(ref_allele.x+ref_allele.y)==0.5));

summary(data0$neighbors.x);
summary(data1$neighbors);
summary(data2$neighbors);
summary(data3$neighbors);

#neighbor_plus <- data3$neighbors+1;
neighbor_plus <- data1_alt$neighbors+1;
neighbor_plus <- data0$neighbors.x+1;
complete_data <- cbind(data1_alt,neighbor_plus);
complete_data <- cbind(data3,neighbor_plus);

#par(mfrow=c(2,2));

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_36b_alt_v1_box.pdf");

#boxplot(ref_allele/(ref_allele+alt_allele) ~ neighbor_plus,data = complete_data, varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,12),ylim=c(0,1),pars = list(outpch=19,cex=0.2));

boxplot(alt_allele/(ref_allele+alt_allele) ~ neighbor_plus,data = complete_data, varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,12),ylim=c(0,1),pars = list(outpch=19,cex=0.2));

abline(h=0.5,lty=2,col="red");

dev.off();

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_36b_v0_box.pdf");

boxplot(ref_allele.x/(ref_allele.x+ref_allele.y) ~ neighbor_plus,data = complete_data, varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,12),ylim=c(0,1),pars = list(outpch=19,cex=0.2));

abline(h=0.5,lty=2,col="red");

dev.off();

###############
# mappability #
###############

library(VennDiagram);

# 1 mismatch

mappability_1mm <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/SNPs_line_40.dm3_ref.l36_m1.mappability.txt",header=TRUE,sep="\t");

data <- merge(data1,mappability_1mm,by.x=c("chr","pos"),by.y=c("chr","position"));

nrow(subset(data,neighbors >= 1));
nrow(subset(data,sum/length != 1));
nrow(subset(data,neighbors >= 1 & sum/length != 1));

nrow(subset(data,neighbors < 1));
nrow(subset(data,sum/length == 1));
nrow(subset(data,neighbors < 1 & sum/length == 1));

area1 <- nrow(subset(data,neighbors >= 1));
area2 <- nrow(subset(data,sum/length != 1));
cross.area <- nrow(subset(data,neighbors >= 1 & sum/length != 1));

venn.plot <- draw.pairwise.venn(area1,area2,cross.area,
	scaled = FALSE,
	fill = c("cyan", "magenta"),
	lty = "blank",
	cex = 2,
	cat.cex = 2,
	cat.col = "black"
);

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/melmel_1mm_venn.pdf");
grid.draw(venn.plot);
dev.off();

nrow(subset(data,ref_allele/(ref_allele+alt_allele) != 0.5 & sum/length != 1))/nrow(subset(data,sum/length != 1))*100;


write.table(data,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_SNPs_1mm.txt",quote=F,sep="\t",row.names=F,col.names=F);

temp <- subset(data,neighbors<1);
write.table(temp,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_SNPs_1mm_unbiased.txt",quote=F,sep="\t",row.names=F,col.names=F);

perfect <- subset(temp,sum/length == 1);
imperfect <- subset(temp,sum/length != 1);

write.table(perfect,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_SNPs_1mm_unbiased_perfect.txt",quote=F,sep="\t",row.names=F,col.names=F);

nrow(subset(temp,ref_allele/(ref_allele+alt_allele) != 0.5 & sum/length != 1))/nrow(subset(temp,ref_allele/(ref_allele+alt_allele) != 0.5))*100; # % unequal ASRA and imperfect mappability 
nrow(subset(temp,ref_allele/(ref_allele+alt_allele) == 0.5 & sum/length != 1))/nrow(subset(temp,ref_allele/(ref_allele+alt_allele) == 0.5))*100; # % equal ASRA and imperfect mappability

# proportion of differentiating sites with perfect mappability and equal allelic abundance
nrow(subset(temp,ref_allele/(ref_allele+alt_allele) != 0.5 & sum/length != 1))/nrow(subset(temp,sum/length != 1))*100;

nonsig <- subset(temp,ref_allele/(ref_allele+alt_allele) == 0.5);
sig <- subset(temp,ref_allele/(ref_allele+alt_allele) != 0.5);
nrow(subset(sig,sum/length < 1))/nrow(sig)*100;
nrow(subset(nonsig,sum/length < 1))/nrow(nonsig)*100;

nonsig_hist <- hist( (nonsig$sum/nonsig$length) ,breaks=seq(0,1,0.05));
sig_hist <- hist( (sig$sum/sig$length) ,breaks=seq(0,1,0.05));
mat <- cbind(nonsig_hist$counts/nrow(nonsig),sig_hist$counts/nrow(sig));
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),xlab="",ylab="",main="",col=c("white","grey"));

perfect_hist <- hist(perfect$ref_allele/(perfect$ref_allele+perfect$alt_allele),breaks=seq(0,1,0.05));
imperfect_hist <- hist(imperfect$ref_allele/(imperfect$ref_allele+imperfect$alt_allele),breaks=seq(0,1,0.05));
mat <- cbind(perfect_hist$counts/nrow(perfect),imperfect_hist$counts/nrow(imperfect));

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_1mm_barplot_by_mapp.pdf");
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),ylim=c(0,1),xlab="",ylab="",main="",col=c("white","grey"));
#par(new=TRUE);
#points(seq(0.05,1,0.05),mat[,1],ylim=c(0,1),col="black",type="b");
#par(new=TRUE);
#points(seq(0.05,1,0.05),mat[,2],ylim=c(0,1),col="grey",type="b");
dev.off();

pdf(file="/Users/kraigrs/Desktop/plots/sim_1mm_pie_by_mapp.pdf");
pie(nrow(perfect),nrow(imperfect),labels="");
dev.off();

pdf(file="/Users/kraigrs/Desktop/plots/sim_1mm_boxplot_perfect.pdf");
boxplot(ref_allele/(ref_allele+alt_allele) ~ neighbor_plus,data = perfect, varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,14),ylim=c(0,1),pars = list(outpch=19,cex=0.4,col=rgb(0,0,0,0.4)));
abline(h = 0.5,col="red",lty="dashed");
dev.off();

pdf(file="/Users/kraigrs/Desktop/plots/sim_1mm_boxplot_imperfect.pdf");
boxplot(ref_allele/(ref_allele+alt_allele) ~ neighbor_plus,data = imperfect, varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,14),ylim=c(0,1),pars = list(outpch=19,cex=0.4,col=rgb(0,0,0,0.4)));
abline(h = 0.5,col="red",lty="dashed");
dev.off();

d1_perfect <- density(perfect$ref_allele/(perfect$ref_allele+perfect$alt_allele),from=0,to=1);
d1_imperfect <- density(imperfect$ref_allele/(imperfect$ref_allele+imperfect$alt_allele),from=0,to=1);

plot(d1_perfect,xlim=c(0,1),ylim=c(0,10),main="",col="blue",xlab="");
par(new=TRUE);
plot(d1_imperfect,xlim=c(0,1),ylim=c(0,10),main="",col="red",xlab="");

# 2 mismatch

mappability_2mm <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/SNPs_line_40.dm3_ref.l36_m2.mappability.txt",header=TRUE,sep="\t");

data <- merge(data2,mappability_2mm,by.x=c("chr","pos"),by.y=c("chr","position"));

nrow(subset(data,neighbors >= 2));
nrow(subset(data,sum/length != 1));
nrow(subset(data,neighbors >= 2 & sum/length != 1));

nrow(data);
nrow(subset(data,neighbors < 2));
nrow(subset(data,sum/length == 1));
nrow(subset(data,neighbors < 2 & sum/length == 1));

nrow(subset(data,ref_allele/(ref_allele+alt_allele) != 0.5 & sum/length != 1))/nrow(subset(data,sum/length != 1))*100;

write.table(data,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_SNPs_2mm.txt",quote=F,sep="\t",row.names=F,col.names=F);

temp <- subset(data,neighbors<2);
write.table(temp,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_SNPs_2mm_unbiased.txt",quote=F,sep="\t",row.names=F,col.names=F);

perfect <- subset(temp,sum/length == 1);
imperfect <- subset(temp,sum/length != 1);

write.table(perfect,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_SNPs_2mm_unbiased_perfect.txt",quote=F,sep="\t",row.names=F,col.names=F);

# proportion of differentiating sites with perfect mappability and equal allelic abundance
nrow(subset(temp,ref_allele/(ref_allele+alt_allele) == 0.5 & sum/length == 1))/nrow(subset(temp,sum/length == 1))*100;

nrow(subset(data,ref_allele/(ref_allele+alt_allele) == 0.5 & sum/length == 1))/nrow(subset(data,ref_allele/(ref_allele+alt_allele) == 0.5))*100; #  
nrow(subset(data,ref_allele/(ref_allele+alt_allele) == 0.5 & sum/length != 1))/nrow(subset(data,ref_allele/(ref_allele+alt_allele) == 0.5))*100; # 

nrow(subset(temp,ref_allele/(ref_allele+alt_allele) != 0.5 & sum/length != 1))/nrow(subset(temp,ref_allele/(ref_allele+alt_allele) != 0.5))*100; # % unequal ASRA and imperfect mappability 
nrow(subset(temp,ref_allele/(ref_allele+alt_allele) == 0.5 & sum/length != 1))/nrow(subset(temp,ref_allele/(ref_allele+alt_allele) == 0.5))*100; # % equal ASRA and imperfect mappability

nonsig <- subset(temp,ref_allele/(ref_allele+alt_allele) == 0.5);
sig <- subset(temp,ref_allele/(ref_allele+alt_allele) != 0.5);
nrow(subset(sig,sum/length < 1))/nrow(sig)*100;
nrow(subset(nonsig,sum/length < 1))/nrow(nonsig)*100;

nonsig_hist <- hist( (nonsig$sum/nonsig$length) ,breaks=seq(0,1,0.05));
sig_hist <- hist( (sig$sum/sig$length) ,breaks=seq(0,1,0.05));
mat <- cbind(nonsig_hist$counts/nrow(nonsig),sig_hist$counts/nrow(sig));
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),xlab="",ylab="",main="",col=c("white","grey"));

perfect_hist <- hist(perfect$ref_allele/(perfect$ref_allele+perfect$alt_allele),breaks=seq(0,1,0.05));
imperfect_hist <- hist(imperfect$ref_allele/(imperfect$ref_allele+imperfect$alt_allele),breaks=seq(0,1,0.05));
mat <- cbind(perfect_hist$counts/nrow(perfect),imperfect_hist$counts/nrow(imperfect));

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_2mm_barplot_by_mapp.pdf");
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),ylim=c(0,1),xlab="",ylab="",main="",col=c("white","grey"));
#par(new=TRUE);
#points(seq(0.05,1,0.05),mat[,1],ylim=c(0,1),col="black",type="b");
#par(new=TRUE);
#points(seq(0.05,1,0.05),mat[,2],ylim=c(0,1),col="grey",type="b");
dev.off();

pdf(file="/Users/kraigrs/Desktop/plots/sim_2mm_boxplot_perfect.pdf");
boxplot(ref_allele/(ref_allele+alt_allele) ~ neighbor_plus,data = perfect, varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,14),ylim=c(0,1),pars = list(outpch=19,cex=0.4,col=rgb(0,0,0,0.4)));
abline(h = 0.5,col="red",lty="dashed");
dev.off();

pdf(file="/Users/kraigrs/Desktop/plots/sim_2mm_boxplot_imperfect.pdf");
boxplot(ref_allele/(ref_allele+alt_allele) ~ neighbor_plus,data = imperfect, varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,14),ylim=c(0,1),pars = list(outpch=19,cex=0.4,col=rgb(0,0,0,0.4)));
abline(h = 0.5,col="red",lty="dashed");
dev.off();

pdf(file="/Users/kraigrs/Desktop/plots/sim_2mm_pie_by_mapp.pdf");
pie(nrow(perfect),nrow(imperfect),labels="");
dev.off();

d2_perfect <- density(perfect$ref_allele/(perfect$ref_allele+perfect$alt_allele));
d2_imperfect <- density(imperfect$ref_allele/(imperfect$ref_allele+imperfect$alt_allele));

plot(d2_perfect,xlim=c(0,1),ylim=c(0,10),main="",col="blue",xlab="");
par(new=TRUE);
plot(d2_imperfect,xlim=c(0,1),ylim=c(0,10),main="",col="red",xlab="");

# 3 mismatch

mappability_3mm <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/SNPs_line_40.dm3_ref.l36_m3.mappability.txt",header=TRUE,sep="\t");

data <- merge(data3,mappability_3mm,by.x=c("chr","pos"),by.y=c("chr","position"));

nrow(data);

nrow(subset(data,neighbors >= 3));
nrow(subset(data,sum/length != 1));
nrow(subset(data,neighbors >= 3 & sum/length != 1));

nrow(subset(data,neighbors < 3));
nrow(subset(data,sum/length == 1));
nrow(subset(data,neighbors < 3 & sum/length == 1));

nrow(subset(data,ref_allele/(ref_allele+alt_allele) != 0.5 & sum/length != 1))/nrow(subset(data,sum/length != 1))*100;

write.table(data,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_SNPs_3mm.txt",quote=F,sep="\t",row.names=F,col.names=F);

temp <- subset(data,neighbors<3);
write.table(temp,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_SNPs_3mm_unbiased.txt",quote=F,sep="\t",row.names=F,col.names=F);

perfect <- subset(temp,sum/length == 1);
imperfect <- subset(temp,sum/length != 1);

write.table(perfect,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_SNPs_3mm_unbiased_perfect.txt",quote=F,sep="\t",row.names=F,col.names=F);

# proportion of differentiating sites with perfect mappability and equal allelic abundance
nrow(subset(temp,ref_allele/(ref_allele+alt_allele) == 0.5 & sum/length == 1))/nrow(subset(temp,sum/length == 1))*100;

nrow(subset(temp,ref_allele/(ref_allele+alt_allele) != 0.5 & sum/length != 1))/nrow(subset(temp,ref_allele/(ref_allele+alt_allele) != 0.5))*100; # % unequal ASRA and imperfect mappability 
nrow(subset(temp,ref_allele/(ref_allele+alt_allele) == 0.5 & sum/length != 1))/nrow(subset(temp,ref_allele/(ref_allele+alt_allele) == 0.5))*100; # % equal ASRA and imperfect mappability

perfect <- subset(temp,sum/length == 1);
imperfect <- subset(temp,sum/length != 1);

nonsig <- subset(temp,ref_allele/(ref_allele+alt_allele) == 0.5);
sig <- subset(temp,ref_allele/(ref_allele+alt_allele) != 0.5);
nrow(subset(sig,sum/length < 1))/nrow(sig)*100;
nrow(subset(nonsig,sum/length < 1))/nrow(nonsig)*100;

nonsig_hist <- hist( (nonsig$sum/nonsig$length) ,breaks=seq(0,1,0.05));
sig_hist <- hist( (sig$sum/sig$length) ,breaks=seq(0,1,0.05));
mat <- cbind(nonsig_hist$counts/nrow(nonsig),sig_hist$counts/nrow(sig));
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),xlab="",ylab="",main="",col=c("white","grey"));

perfect_hist <- hist(perfect$ref_allele/(perfect$ref_allele+perfect$alt_allele),breaks=seq(0,1,0.05));
imperfect_hist <- hist(imperfect$ref_allele/(imperfect$ref_allele+imperfect$alt_allele),breaks=seq(0,1,0.05));
mat <- cbind(perfect_hist$counts/nrow(perfect),imperfect_hist$counts/nrow(imperfect));

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_3mm_barplot_by_mapp.pdf");
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),ylim=c(0,1),xlab="",ylab="",main="",col=c("white","grey"));
#par(new=TRUE);
#points(seq(0.05,1,0.05),mat[,1],ylim=c(0,1),col="black",type="b");
#par(new=TRUE);
#points(seq(0.05,1,0.05),mat[,2],ylim=c(0,1),col="grey",type="b");
dev.off();

pdf(file="/Users/kraigrs/Desktop/plots/sim_3mm_pie_by_mapp.pdf");
pie(nrow(perfect),nrow(imperfect),labels="");
dev.off();

pdf(file="/Users/kraigrs/Desktop/plots/sim_3mm_boxplot_perfect.pdf");
boxplot(ref_allele/(ref_allele+alt_allele) ~ neighbor_plus,data = perfect, varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,14),ylim=c(0,1),pars = list(outpch=19,cex=0.4,col=rgb(0,0,0,0.4)));
abline(h = 0.5,col="red",lty="dashed");
dev.off();

pdf(file="/Users/kraigrs/Desktop/plots/sim_3mm_boxplot_imperfect.pdf");
boxplot(ref_allele/(ref_allele+alt_allele) ~ neighbor_plus,data = imperfect, varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,14),ylim=c(0,1),pars = list(outpch=19,cex=0.4,col=rgb(0,0,0,0.4)));
abline(h = 0.5,col="red",lty="dashed");
dev.off();

d3_perfect <- density(perfect$ref_allele/(perfect$ref_allele+perfect$alt_allele));
d3_imperfect <- density(imperfect$ref_allele/(imperfect$ref_allele+imperfect$alt_allele));

plot(d3_perfect,xlim=c(0,1),ylim=c(0,10),main="",col="blue",xlab="");
par(new=TRUE);
plot(d3_imperfect,xlim=c(0,1),ylim=c(0,10),main="",col="red",xlab="");

# 0 mismatches

ref_mappability <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/SNPs_line_40.dm3_ref.l36_m0.mappability.txt",header=TRUE,sep="\t");
alt_mappability <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/SNPs_line_40.dm3_alt_line_40.l36_m0.mappability.txt",header=TRUE,sep="\t");
mappability_0mm <- merge(ref_mappability,alt_mappability,by.x=c("chr","position"),by.y=c("chr","position"));

data <- merge(data0,mappability_0mm,by.x=c("chr","pos"),by.y=c("chr","position"));

nrow(data);
nrow(subset(data,((sum.x/length.x)+(sum.y/length.y)) == 2));


nrow(subset(data,ref_allele.x/(ref_allele.x+ref_allele.y) != 0.5 & ((sum.x/length.x)+(sum.y/length.y)) != 2))/nrow(subset(data,((sum.x/length.x)+(sum.y/length.y)) != 2))*100;


write.table(data,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_SNPs_0mm.txt",quote=F,sep="\t",row.names=F,col.names=F);

write.table(data,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_SNPs_0mm_unbiased.txt",quote=F,sep="\t",row.names=F,col.names=F);

perfect <- subset(data,((sum.x/length.x)+(sum.y/length.y)) == 2);
imperfect <- subset(data,((sum.x/length.x)+(sum.y/length.y)) != 2);

write.table(perfect,file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_SNPs_0mm_unbiased_perfect.txt",quote=F,sep="\t",row.names=F,col.names=F);

# proportion of differentiating sites with perfect mappability and equal allelic abundance
nrow(subset(data,ref_allele.x/(ref_allele.x+ref_allele.y) == 0.5 & ((sum.x/length.x)+(sum.y/length.y)) == 2))/nrow(subset(data,((sum.x/length.x)+(sum.y/length.y)) == 2))*100;



perfect_hist <- hist(perfect$ref_allele.x/(perfect$ref_allele.x+perfect$ref_allele.y),breaks=seq(0,1,0.05));
imperfect_hist <- hist(imperfect$ref_allele.x/(imperfect$ref_allele.x+imperfect$ref_allele.y),breaks=seq(0,1,0.05));
mat <- cbind(perfect_hist$counts/nrow(perfect),imperfect_hist$counts/nrow(imperfect));

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/sim_0mm_barplot_by_mapp.pdf");
barplot(t(mat),beside=TRUE,names.arg=seq(0.05,1,0.05),ylim=c(0,1),xlab="",ylab="",main="",col=c("white","grey"));
#par(new=TRUE);
#points(seq(0.05,1,0.05),mat[,1],ylim=c(0,1),col="black",type="b");
#par(new=TRUE);
#points(seq(0.05,1,0.05),mat[,2],ylim=c(0,1),col="grey",type="b");
dev.off();

pdf(file="/Users/kraigrs/Desktop/plots/sim_0mm_pie_by_mapp.pdf");
pie(nrow(perfect),nrow(imperfect),labels="");
dev.off();

pdf(file="/Users/kraigrs/Desktop/plots/sim_0mm_boxplot_perfect.pdf");
boxplot(dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele) ~ neighbor_plus,data = perfect, varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,14),ylim=c(0,1),pars = list(outpch=19,cex=0.4,col=rgb(0,0,0,0.4)));
abline(h = 0.5,col="red",lty="dashed");
dev.off();

pdf(file="/Users/kraigrs/Desktop/plots/sim_0mm_boxplot_imperfect.pdf");
boxplot(dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele) ~ neighbor_plus,data = imperfect, varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,14),ylim=c(0,1),pars = list(outpch=19,cex=0.4,col=rgb(0,0,0,0.4)));
abline(h = 0.5,col="red",lty="dashed");
dev.off();

d0_perfect <- density(perfect$dm3_ref_ref_allele/(perfect$dm3_ref_ref_allele+perfect$dm3_alt_alt_allele));
d0_imperfect <- density(imperfect$dm3_ref_ref_allele/(imperfect$dm3_ref_ref_allele+imperfect$dm3_alt_alt_allele));

plot(d0_perfect,xlim=c(0,1),ylim=c(0,10),main="",col="blue",xlab="");
par(new=TRUE);
plot(d0_imperfect,xlim=c(0,1),ylim=c(0,10),main="",col="red",xlab="");


#########################

plot((data$sum.x/data$length.x)/(data$sum.x/data$length.x+data$sum.y/data$length.y),data$dm3_ref_ref_allele/(data$dm3_ref_ref_allele+data$dm3_alt_alt_allele),xlab="mappability",ylab="proportion of reference allele",pch=19,col=rgb(0,0,0,0.5),cex=0.5)

data <- subset(data,dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele) != 0.5 & (sum.x/length.x)/(sum.x/length.x+sum.y/length.y) != 0.5);

data <- subset(data,dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele) != 0.5);

nrow(subset(data, sign((sum.x/length.x)/(sum.x/length.x+sum.y/length.y)-0.5) == sign(dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele)-0.5) ));

multiple <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.bowtie_v0_m1.SNPs.txt",header=TRUE,sep="\t");

data <- merge(multiple,SNP_mappability,by.x=c("chr","pos"),by.y=c("chr","position"));

nrow(subset(data,dm3_ref_ref_allele > 0 & dm3_alt_alt_allele > 0 & log2(dm3_ref_ref_allele/dm3_alt_alt_allele) != 0));
nrow(subset(data,dm3_ref_ref_allele > 0 & dm3_alt_alt_allele > 0 & log2(dm3_ref_ref_allele/dm3_alt_alt_allele) != 0 & log2((sum.x/length.x)/(sum.y/length.y)) != 0));

nrow(subset(data,dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele) != 0.5));
nrow(subset(data,dm3_ref_ref_allele/(dm3_ref_ref_allele+dm3_alt_alt_allele) != 0.5 & log2((sum.x/length.x)/(sum.y/length.y)) != 0));

poor_map <- subset(data,dm3_ref_ref_allele > 0 & dm3_alt_alt_allele > 0 & log2(dm3_ref_ref_allele/dm3_alt_alt_allele) != 0 & log2((sum.x/length.x)/(sum.y/length.y)) != 0);
nrow(subset(poor_map,sign(log2(dm3_ref_ref_allele/dm3_alt_alt_allele)) == sign(log2((sum.x/length.x)/(sum.y/length.y)))));

AI <- subset(data,ref_allele.x > 0 & alt_allele.x > 0 & log2(ref_allele.x/alt_allele.x) != 0);
AI_prop <- nrow(subset(AI,sum/length != 1))/nrow(AI);

ASE <- subset(data,ref_allele.x > 0 & alt_allele.x > 0 & log2(ref_allele.x/alt_allele.x) == 0);
ASE_prop <- nrow(subset(ASE,sum/length != 1))/nrow(ASE);

# SNPs with no neighbors and AI

temp1 <- subset(AI,neighbor==0);
AI_prop <- nrow(subset(temp1,sum/length != 1))/nrow(temp1); 

temp2 <- subset(ASE,neighbor==0);
ASE_prop <- nrow(subset(temp2,sum/length != 1))/nrow(temp2); 

props <- c(AI_prop,ASE_prop);
barplot(props,names.arg=c("AI\nn = 179","not AI\nn = 38,819"),ylim=c(0,1),ylab="Proportion filtered out in each category",main="Filtering imperfect mappability retains many SNPs");

barplot(c(ASE_prop,AI_prop),xlim=c(0,1),horiz=TRUE,names.arg=c("No AI\nn = 179","AI\nn = 38,819"),width=0.3,xlab="Proportion removed due to imperfect mappability",main="SNPs");

##################################
# plot mappability across genome #
##################################

library(lattice);

dm3_ref_exon_mappability <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons.dm3_ref.l50_m0.mappability.txt",header=TRUE,sep="\t");

dm3_alt_exon_mappability <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons.dm3_alt_line_40.l50_m0.mappability.txt",header=TRUE,sep="\t");

exon_mappability = merge(dm3_ref_exon_mappability,dm3_alt_exon_mappability,by.x="locus",by.y="locus");
dm3_ref_avg <- exon_mappability$sum.x/exon_mappability$length.x;
dm3_alt_avg <- exon_mappability$sum.y/exon_mappability$length.y;

mappability <- cbind(exon_mappability,dm3_ref_avg,dm3_alt_avg);

multiple <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.bowtie_v0_m1.exons.txt",header=TRUE,sep="\t");

SNPsInExons <- read.table("/Users/kraigrs/Wittkopp/DGRP/DGRP_line_40_SNPs_in_const.txt",header=FALSE,sep="\t");

temp1 <- merge(multiple,mappability,by.x="gene_exon",by.y="locus");
temp2 <- merge(temp1,SNPsInExons,by.x="gene_exon",by.y="V8");

temp3 <- subset(temp2,dm3 > 0 & line_40 > 0);

xyplot(dm3_ref_avg/dm3_alt_avg ~ (V2+V3)/2 | V1, data = temp3,xlab="Midpoint of exon",ylab="dm3/line_40 mappability",main="Mappability across both allele-specific genomes",pch=19,cex=0.4,col=rgb(0,0,0,0.3),layout=c(2,3));

xyplot(dm3_ref_avg/dm3_alt_avg ~ log2(dm3/line_40) | V1, data = temp3,xlab="log2(dm3/line_40) ASE",ylab="dm3/line_40 mappability",main="Concordance of mappability direction and bias in ASE",pch=19,cex=0.4,col=rgb(0,0,0,0.3));

# how many exons show opposite signs of mappability and bias? 0.2%, so 99.8% show same sign, awesome!
# also, of the 660 that do show differential mappability, 237 favor dm3_ref and 423 favor dm3_alt
sum( sign(log2(temp3$dm3_ref_avg/temp3$dm3_alt_avg)) != sign(log2(temp3$dm3/temp3$line_40)) );

# correlation? r^2 = 0.6
cor(log2(temp3$dm3_ref_avg/temp3$dm3_alt_avg),log2(temp3$dm3_ref/temp3$dm3_alt))^2;

###############

SNPsInExons <- read.table("/Users/kraigrs/Wittkopp/DGRP/DGRP_line_40_SNPs_in_const.txt",header=FALSE,sep="\t");

v0 <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.bowtie_v0_m1.exons.txt",header=TRUE,sep="\t");
v1 <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v1_m1.exons.txt",header=TRUE,sep="\t");
v2 <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v2_m1.exons.txt",header=TRUE,sep="\t");
v3 <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v3_m1.exons.txt",header=TRUE,sep="\t");

v0 <- merge(v0,SNPsInExons,by.x="gene_exon",by.y="V8");
v1 <- merge(v1,SNPsInExons,by.x="gene_exon",by.y="V8");
v2 <- merge(v2,SNPsInExons,by.x="gene_exon",by.y="V8");
v3 <- merge(v3,SNPsInExons,by.x="gene_exon",by.y="V8");

nrow(subset(v0,dm3>0&line_40>0&log2(dm3/line_40)!=0&V7>0))/nrow(v0);
nrow(subset(v1,dm3>0&line_40>0&log2(dm3/line_40)!=0&V7>0))/nrow(v1);
nrow(subset(v2,dm3>0&line_40>0&log2(dm3/line_40)!=0&V7>0))/nrow(v2);
nrow(subset(v3,dm3>0&line_40>0&log2(dm3/line_40)!=0&V7>0))/nrow(v3);

nrow(subset(v0,dm3>0&line_40>0&V7>0));
nrow(subset(v1,dm3>0&line_40>0&V7>0));
nrow(subset(v2,dm3>0&line_40>0&V7>0));
nrow(subset(v3,dm3>0&line_40>0&V7>0));


# mappability for single-genome approach

data <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_ref.bowtie_v1_m1.SNPs.txt",header=TRUE,sep="\t");

SNP_mappability <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/SNPs_line_40.dm3_ref.l50_m1.mappability.txt",header=TRUE,sep="\t");

data <- merge(complete_data,SNP_mappability,by.x=c("chr","pos"),by.y=c("chr","position"));

AI <- subset(data, ref_allele/(ref_allele+alt_allele) != 0.5));
AI_prop <- nrow(subset(AI,sum/length != 1))/nrow(AI);

ASE <- subset(data,ref_allele.x > 0 & alt_allele.x > 0 & log2(ref_allele.x/alt_allele.x) == 0);
ASE_prop <- nrow(subset(ASE,sum/length != 1))/nrow(ASE);

# SNPs with no neighbors and AI

temp1 <- subset(AI,neighbor==0);
AI_prop <- nrow(subset(temp1,sum/length != 1))/nrow(temp1); 

temp2 <- subset(ASE,neighbor==0);
ASE_prop <- nrow(subset(temp2,sum/length != 1))/nrow(temp2);

##############
# alternative allele used as genome

data <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/line_40/constExons_single_bp50_error0_tiled_line_40.dm3_alt_line_40.bowtie_v1_m1.SNPs.txt",header=TRUE,sep="\t");
neighbors_plus <- data$neighbors+1;
data <- cbind(data,neighbors_plus);

boxplot(alt_allele/(alt_allele+ref_allele) ~ neighbors_plus,data = data, varwidth = TRUE,xlab="",ylab="",main="",xlim=c(1,14),ylim=c(0,1),pars = list(outpch=19,cex=0.2));

abline(h=0.5,lty=2,col="red");

pie(c(nrow(subset(data,alt_allele/(ref_allele+alt_allele)!=0.5)),nrow(data)),col=c("gray","white"),labels="");
