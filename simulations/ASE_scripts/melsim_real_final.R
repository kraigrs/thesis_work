###### mappability ######

data1_tests <- read.table("/Users/kraigrs/Wittkopp/Graze/SRR389095.mate1.berlin-updated-exonic-regions.bowtie_v1_m1.SNPs.binom_tests.txt",header=TRUE,sep="\t");

data2_tests <- read.table("/Users/kraigrs/Wittkopp/Graze/SRR389095.mate1.berlin-updated-exonic-regions.bowtie_v2_m1.SNPs.binom_tests.txt",header=TRUE,sep="\t");

data3_tests <- read.table("/Users/kraigrs/Wittkopp/Graze/SRR389095.mate1.berlin-updated-exonic-regions.bowtie_v3_m1.SNPs.binom_tests.txt",header=TRUE,sep="\t");

data0_tests <- read.table("/Users/kraigrs/Wittkopp/Graze/SRR389095.mate1.berlin_c1674_fsa_masked.bowtie_v0_m1.SNPs.binom_tests.txt",header=TRUE,sep="\t");

map1 <- read.table("/Users/kraigrs/Wittkopp/Graze/berlin-updated-exonic-regions.SNPs.l36_m1.mappability.txt",header=TRUE,sep="\t");
map2 <- read.table("/Users/kraigrs/Wittkopp/Graze/berlin-updated-exonic-regions.SNPs.l36_m2.mappability.txt",header=TRUE,sep="\t");
map3 <- read.table("/Users/kraigrs/Wittkopp/Graze/berlin-updated-exonic-regions.SNPs.l36_m3.mappability.txt",header=TRUE,sep="\t");

berlin_map <- read.table("/Users/kraigrs/Wittkopp/Graze/berlin-updated-exonic-regions_fsa_masked.SNPs.mappability.txt",header=TRUE,sep="\t");
c1674_map <- read.table("/Users/kraigrs/Wittkopp/Graze/c1674-updated-exonic-regions_fsa_masked.SNPs.mappability.txt",header=TRUE,sep="\t");
map0 <- merge(berlin_map,c1674_map,by.x=c("locus","position"),by.y=c("locus","position"));

temp <- read.table("/Users/kraigrs/Wittkopp/Graze/berlin_c1674_fsa.indels.txt",header=TRUE,sep="\t");
key <- read.table("/Users/kraigrs/Wittkopp/Graze/key",header=TRUE,sep="\t");
indels <- merge(temp,key,by.x=c("chr","pos"),by.y=c("locus","fsa"));

relevant <- read.table("/Users/kraigrs/Wittkopp/Graze/berlin-updated-exonic-regions_fsa.36b_windows_SNPs.neighbors.txt",header=TRUE,sep="\t");

relevant <- merge(relevant,key,by.x=c("chr","pos"),by.y=c("locus","mel"));
summary(relevant$relevant);

data1_tests <- merge(data1_tests,relevant,by.x=c("chr","pos"),by.y=c("chr","pos"));
data2_tests <- merge(data2_tests,relevant,by.x=c("chr","pos"),by.y=c("chr","pos"));
data3_tests <- merge(data3_tests,relevant,by.x=c("chr","pos"),by.y=c("chr","pos"));


data1_tests_map1 <- merge(data1_tests,map1,by.x=c("chr","pos"),by.y=c("locus","position"));
data2_tests_map2 <- merge(data2_tests,map2,by.x=c("chr","pos"),by.y=c("locus","position"));
data3_tests_map3 <- merge(data3_tests,map3,by.x=c("chr","pos"),by.y=c("locus","position"));
data0_tests_map0 <- merge(data0_tests,map0,by.x=c("chr","pos"),by.y=c("locus","position"));

data1_tests_map1_indels <- merge(data1_tests_map1,indels,by.x=c("chr","pos"),by.y=c("chr","mel"));
data1_tests_map1_indels <- data1_tests_map1_indels[,c(1:14,16:18)];

data2_tests_map2_indels <- merge(data2_tests_map2,indels,by.x=c("chr","pos"),by.y=c("chr","mel"));
data2_tests_map2_indels <- data2_tests_map2_indels[,c(1:14,16:18)];

data3_tests_map3_indels <- merge(data3_tests_map3,indels,by.x=c("chr","pos"),by.y=c("chr","mel"));
data3_tests_map3_indels <- data3_tests_map3_indels[,c(1:14,16:18)];

data0_tests_map0_indels <- merge(data0_tests_map0,indels,by.x=c("chr","pos"),by.y=c("chr","pos"));


# 1 mismatch

nrow(data1_tests_map1_indels);
nrow(subset(data1_tests_map1_indels,qvals >= 0.05));

nrow(subset(data1_tests_map1_indels,relevant <= 1));
nrow(subset(data1_tests_map1_indels,relevant <= 1 & qvals >= 0.05));

nrow(subset(data1_tests_map1_indels,relevant <= 1 & sum/length == 1));
nrow(subset(data1_tests_map1_indels,relevant <= 1 & sum/length == 1 & qvals >= 0.05));

nrow(subset(data1_tests_map1_indels,relevant <= 1 & sum/length == 1 & ref_indel + alt_indel == 0));
nrow(subset(data1_tests_map1_indels,relevant <= 1 & sum/length == 1 & ref_indel + alt_indel == 0 & qvals >= 0.05));

# 2 mismatches

nrow(data2_tests_map2_indels);
nrow(subset(data2_tests_map2_indels,qvals >= 0.05));

nrow(subset(data2_tests_map2_indels,relevant <= 2));
nrow(subset(data2_tests_map2_indels,relevant <= 2 & qvals >= 0.05));

nrow(subset(data2_tests_map2_indels,relevant <= 2 & sum/length == 1));
nrow(subset(data2_tests_map2_indels,relevant <= 2 & sum/length == 1 & qvals >= 0.05));

nrow(subset(data2_tests_map2_indels,relevant <= 2 & sum/length == 1 & ref_indel + alt_indel == 0));
nrow(subset(data2_tests_map2_indels,relevant <= 2 & sum/length == 1 & ref_indel + alt_indel == 0 & qvals >= 0.05));


# 3 mismatches

nrow(data3_tests_map3_indels);
nrow(subset(data3_tests_map3_indels,qvals >= 0.05));

nrow(subset(data3_tests_map3_indels,relevant <= 3));
nrow(subset(data3_tests_map3_indels,relevant <= 3 & qvals >= 0.05));

nrow(subset(data3_tests_map3_indels,relevant <= 3 & sum/length == 1));
nrow(subset(data3_tests_map3_indels,relevant <= 3 & sum/length == 1 & qvals >= 0.05));

nrow(subset(data3_tests_map3_indels,relevant <= 3 & sum/length == 1 & ref_indel + alt_indel == 0));
nrow(subset(data3_tests_map3_indels,relevant <= 3 & sum/length == 1 & ref_indel + alt_indel == 0 & qvals >= 0.05));

# 0 mismatches

nrow(data0_tests_map0_indels);
nrow(subset(data0_tests_map0_indels,qvals >= 0.05));

nrow(subset(data0_tests_map0_indels,((sum.x/length.x)+(sum.y/length.y)) == 2));
nrow(subset(data0_tests_map0_indels,((sum.x/length.x)+(sum.y/length.y)) == 2 & qvals >= 0.05));

nrow(subset(data0_tests_map0_indels,((sum.x/length.x)+(sum.y/length.y)) == 2 & ref_indel + alt_indel == 0));
nrow(subset(data0_tests_map0_indels,((sum.x/length.x)+(sum.y/length.y)) == 2 & ref_indel + alt_indel == 0 & qvals >= 0.05));


temp1 <- subset(data1_tests_map1_indels,relevant <= 1 & sum/length == 1 & ref_indel + alt_indel == 0);
temp2 <- subset(data2_tests_map2_indels,relevant <= 2 & sum/length == 1 & ref_indel + alt_indel == 0);
temp3 <- subset(data3_tests_map3_indels,relevant <= 3 & sum/length == 1 & ref_indel + alt_indel == 0);
#temp3 <- data3_tests;
temp0 <- subset(data0_tests_map0,((sum.x/length.x)+(sum.y/length.y)) == 2);

# 1 mismatch vs. parental genomes

compare <- merge(temp0,temp1,by.x=c("chr","pos"),by.y=c("chr","fsa"));

cor(compare$berlin_ref_allele/(compare$berlin_ref_allele+compare$c1674_alt_allele) , compare$ref_allele/(compare$ref_allele+compare$alt_allele))^2;

nrow(compare);
s_ns <- subset(compare,qvals.x < 0.05 & qvals.y >= 0.05);
ns_s <- subset(compare,qvals.x >= 0.05 & qvals.y < 0.05);
tot <- nrow(s_ns)+nrow(ns_s);
a <- nrow(subset(s_ns,ref_allele/(ref_allele+alt_allele) > berlin_ref_allele/(berlin_ref_allele+c1674_alt_allele)));
b <- nrow(subset(ns_s,ref_allele/(ref_allele+alt_allele) > berlin_ref_allele/(berlin_ref_allele+c1674_alt_allele)));
(a+b)/tot;

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/compare_real_data_1mm.pdf");
plot( compare$berlin_ref_allele/(compare$berlin_ref_allele+compare$c1674_alt_allele) , compare$ref_allele/(compare$ref_allele+compare$alt_allele) ,xlim=c(0,1),ylim=c(0,1),pch=19,cex=0.5,col=rgb(0,0,0,0.05),xlab="",ylab="",main="");
dev.off();

n_n <- subset(compare,qvals.x >= 0.05 & qvals.y >= 0.05);
s_ns <- subset(compare,qvals.x < 0.05 & qvals.y >= 0.05);
ns_s <- subset(compare,qvals.x >= 0.05 & qvals.y < 0.05);
s_s <- subset(compare,qvals.x < 0.05 & qvals.y < 0.05);

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/pie_real_data_1mm.pdf");
pie(c(nrow(n_n),nrow(s_ns),nrow(ns_s),nrow(s_s)),col=c("grey","red","blue","purple"),labels="");
dev.off();

# 2 mismatches vs. parental genomes

compare <- merge(temp0,temp2,by.x=c("chr","pos"),by.y=c("chr","fsa"));

cor(compare$berlin_ref_allele/(compare$berlin_ref_allele+compare$c1674_alt_allele) , compare$ref_allele/(compare$ref_allele+compare$alt_allele))^2;

nrow(compare);
s_ns <- subset(compare,qvals.x < 0.05 & qvals.y >= 0.05);
ns_s <- subset(compare,qvals.x >= 0.05 & qvals.y < 0.05);
tot <- nrow(s_ns)+nrow(ns_s);
a <- nrow(subset(s_ns,ref_allele/(ref_allele+alt_allele) > berlin_ref_allele/(berlin_ref_allele+c1674_alt_allele)));
b <- nrow(subset(ns_s,ref_allele/(ref_allele+alt_allele) > berlin_ref_allele/(berlin_ref_allele+c1674_alt_allele)));
(a+b)/tot;

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/compare_real_data_2mm.pdf");
plot( compare$berlin_ref_allele/(compare$berlin_ref_allele+compare$c1674_alt_allele) , compare$ref_allele/(compare$ref_allele+compare$alt_allele) ,xlim=c(0,1),ylim=c(0,1),pch=19,cex=0.5,col=rgb(0,0,0,0.05),xlab="",ylab="",main="");
dev.off();

n_n <- subset(compare,qvals.x >= 0.05 & qvals.y >= 0.05);
s_ns <- subset(compare,qvals.x < 0.05 & qvals.y >= 0.05);
ns_s <- subset(compare,qvals.x >= 0.05 & qvals.y < 0.05);
s_s <- subset(compare,qvals.x < 0.05 & qvals.y < 0.05);

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/pie_real_data_2mm.pdf");
pie(c(nrow(n_n),nrow(s_ns),nrow(ns_s),nrow(s_s)),col=c("grey","red","blue","purple"),labels="");
dev.off();


# 3 mismatches vs. parental genomes

compare <- merge(temp0,temp3,by.x=c("chr","pos"),by.y=c("chr","fsa"));

cor(compare$berlin_ref_allele/(compare$berlin_ref_allele+compare$c1674_alt_allele) , compare$ref_allele/(compare$ref_allele+compare$alt_allele))^2;

nrow(compare);
s_ns <- subset(compare,qvals.x < 0.05 & qvals.y >= 0.05);
ns_s <- subset(compare,qvals.x >= 0.05 & qvals.y < 0.05);
tot <- nrow(s_ns)+nrow(ns_s);
a <- nrow(subset(s_ns,ref_allele/(ref_allele+alt_allele) > berlin_ref_allele/(berlin_ref_allele+c1674_alt_allele)));
b <- nrow(subset(ns_s,ref_allele/(ref_allele+alt_allele) > berlin_ref_allele/(berlin_ref_allele+c1674_alt_allele)));
(a+b)/tot;

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/compare_real_data_3mm.pdf");
plot( compare$berlin_ref_allele/(compare$berlin_ref_allele+compare$c1674_alt_allele) , compare$ref_allele/(compare$ref_allele+compare$alt_allele) ,xlim=c(0,1),ylim=c(0,1),pch=19,cex=0.5,col=rgb(0,0,0,0.05),xlab="",ylab="",main="");
dev.off();

n_n <- subset(compare,qvals.x >= 0.05 & qvals.y >= 0.05);
s_ns <- subset(compare,qvals.x < 0.05 & qvals.y >= 0.05);
ns_s <- subset(compare,qvals.x >= 0.05 & qvals.y < 0.05);
s_s <- subset(compare,qvals.x < 0.05 & qvals.y < 0.05);

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/pie_real_data_3mm.pdf");
pie(c(nrow(n_n),nrow(s_ns),nrow(ns_s),nrow(s_s)),col=c("grey","red","blue","purple"),labels="");
dev.off();


# 3 mismatches (biased) vs. parental genomes

compare <- merge(temp0,data3_tests,by.x=c("chr","pos"),by.y=c("chr","fsa"));

cor(compare$berlin_ref_allele/(compare$berlin_ref_allele+compare$c1674_alt_allele) , compare$ref_allele/(compare$ref_allele+compare$alt_allele))^2;

nrow(compare);
s_ns <- subset(compare,qvals.x < 0.05 & qvals.y >= 0.05);
ns_s <- subset(compare,qvals.x >= 0.05 & qvals.y < 0.05);
tot <- nrow(s_ns)+nrow(ns_s);
a <- nrow(subset(s_ns,ref_allele/(ref_allele+alt_allele) > berlin_ref_allele/(berlin_ref_allele+c1674_alt_allele)));
b <- nrow(subset(ns_s,ref_allele/(ref_allele+alt_allele) > berlin_ref_allele/(berlin_ref_allele+c1674_alt_allele)));
(a+b)/tot;

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/compare_real_data_3mm_biased.pdf");
plot( compare$berlin_ref_allele/(compare$berlin_ref_allele+compare$c1674_alt_allele) , compare$ref_allele/(compare$ref_allele+compare$alt_allele) ,xlim=c(0,1),ylim=c(0,1),pch=19,cex=0.5,col=rgb(0,0,0,0.05),xlab="",ylab="",main="");
dev.off();

n_n <- subset(compare,qvals.x >= 0.05 & qvals.y >= 0.05);
s_ns <- subset(compare,qvals.x < 0.05 & qvals.y >= 0.05);
ns_s <- subset(compare,qvals.x >= 0.05 & qvals.y < 0.05);
s_s <- subset(compare,qvals.x < 0.05 & qvals.y < 0.05);

pdf(file="/Users/kraigrs/Wittkopp/LabResearch/Simulating_ASE/revisions/new_plots/pie_real_data_3mm_biased.pdf");
pie(c(nrow(n_n),nrow(s_ns),nrow(ns_s),nrow(s_s)),col=c("grey","red","blue","purple"),labels="");
dev.off();