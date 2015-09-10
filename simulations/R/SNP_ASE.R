#data1 <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/constExons_single_bp50_error0_tiled.dm3_ref.bowtie_mm1.SNP_ASE.txt",header=TRUE,sep="\t");

#data2 <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/constExons_single_bp50_error0_tiled.dm3_alt.bowtie_mm1.SNP_ASE.txt",header=TRUE,sep="\t");

#data <- merge(data1,data2,by.x=c("chr","pos"),by.y=c("chr","pos"));

data1 <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/constExons_single_bp50_error0_tiled.dm3_ref.bowtie_v1_m1.SNP_ASE.txt",header=TRUE,sep="\t");

data2 <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/constExons_single_bp50_error0_tiled.dm3_ref.bowtie_v2_m1.SNP_ASE.txt",header=TRUE,sep="\t");

data3 <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/constExons_single_bp50_error0_tiled.dm3_ref.bowtie_v3_m1.SNP_ASE.txt",header=TRUE,sep="\t");

temp1 <- merge(data1,data2,by.x=c("chr","pos"),by.y=c("chr","pos"));
data <- merge(temp1,data3,by.x=c("chr","pos"),by.y=c("chr","pos"));

complete_data <- NULL;

###### figure out how many SNPs are within 50b of each SNP ######
### start iterative run here

chrom <- subset(data,data$chr == "chrX");
chrom_sorted <- chrom[order(chrom$pos),];
chrom <- chrom_sorted;

steps_right = 0;
steps_left = 0;

SNPs_left <- mat.or.vec(nrow(chrom),1);
SNPs_right <- mat.or.vec(nrow(chrom),1);

for(i in 1:nrow(chrom))
{
	if(i == 1)
	{
		SNPs_left[i] <- 0;
		k <- i + 1;
		right <- chrom$pos[k] - chrom$pos[i];
		
		while(right < 50)
		{
			steps_right <- steps_right + 1;
			if(k == nrow(chrom)){break;}
			k <- k + 1;
			right <- chrom$pos[k] - chrom$pos[i];
		}
		SNPs_right[i] = steps_right;
	}
	
	if(i == nrow(chrom))
	{
		SNPs_right[i] <- 0;
		j <- i - 1;
		left <- chrom$pos[i] - chrom$pos[j];
		
		while(left < 50)
		{
			steps_left <- steps_left + 1;
			if(j == 1){break;}
			j <- j - 1;
			left <- chrom$pos[i] - chrom$pos[j];
		}
		SNPs_left[i] = steps_left;
	}
	
	j <- i - 1;
	k <- i + 1;
	
	steps_left <- 0;
	steps_right <- 0;
	
	if(j != 0 & k != (nrow(chrom)+1))
	{
		left <- chrom$pos[i] - chrom$pos[j];
		right <- chrom$pos[k] - chrom$pos[i];
	
		while(left < 50)
		{
			steps_left <- steps_left + 1;
			j <- j - 1;
			if(j == 0){break;}
			left <- chrom$pos[i] - chrom$pos[j];
		}
		SNPs_left[i] = steps_left;
	
		while(right < 50)
		{
			steps_right <- steps_right + 1;
			k <- k + 1;
			if(k == (nrow(chrom)+1)){break;}
			right <- chrom$pos[k] - chrom$pos[i];
		}
		SNPs_right[i] = steps_right;
	}
}

neighbor <- SNPs_left+SNPs_right;
chrom <- cbind(chrom,neighbor);

complete_data <- rbind(complete_data,chrom);

# run the previous iteratively to get complete dataset

#write.table(complete_data,file="/Users/kraigrs/Wittkopp/Simulations/tiled/constExons_single_bp50_error0_tiled.dm3_ref.bowtie_mm1.neighbor_SNPs.txt",quote=FALSE,sep="\t",row.names=FALSE);

#write.table(ref_bad_SNPs,file="/Users/kraigrs/Wittkopp/Simulations/tiled/ref_bad_SNPs.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE);

#write.table(alt_bad_SNPs,file="/Users/kraigrs/Wittkopp/Simulations/tiled/alt_bad_SNPs.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE);

# 1 mismatch
boxplot(log2(ref_allele.x/alt_allele.x) ~ neighbor,data = subset(complete_data,ref_allele.x > 0 & alt_allele.x > 0),varwidth = TRUE,xlab="Number of neighboring SNPs",ylab="Distribution of log2(ref/alt)",border=rgb(0,0,0,0.5),xlim=c(0,13),ylim=c(-6,6));

par(new=TRUE);

# 2 mismatches
boxplot(log2(ref_allele.y/alt_allele.y) ~ neighbor,data = subset(complete_data, ref_allele.y > 0 & alt_allele.y > 0),varwidth = TRUE,xlab="Number of neighboring SNPs",ylab="Distribution of log2(ref/alt)",border=rgb(1,0,0,0.5),xlim=c(0,13),ylim=c(-6,6));

par(new=TRUE);

# 3 mismatches
boxplot(log2(ref_allele/alt_allele) ~ neighbor,data = subset(complete_data, ref_allele > 0 & alt_allele > 0),varwidth = TRUE,xlab="Number of neighboring SNPs",ylab="Distribution of log2(ref/alt)",border=rgb(0,0,1,0.5),xlim=c(0,13),ylim=c(-6,6));

#############################################
# plot medians for the different mismatches #
#############################################

summary(complete_data$neighbor); # range of neighboring SNPs [0,22]

multiple <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/constExons_single_bp50_error0_tiled.bowtie.SNP_ASE.txt",header=TRUE,sep="\t");

both <- merge(complete_data,multiple,by.x=c("chr","pos"),by.y=c("chr","pos"));

medians <- mat.or.vec(14,5);
colnames(medians) <- c("neighbor","v1","v2","v3","multiple");

vars <- mat.or.vec(14,5);
colnames(vars) <- c("neighbor","v1","v2","v3","multiple");

for(i in 1:14)
{
	j <- i-1;
	temp1 <- subset(both,ref_allele.x > 0 & alt_allele.x > 0 & neighbor == j);
	temp2 <- subset(both,ref_allele.y > 0 & alt_allele.y > 0 & neighbor == j);
	temp3 <- subset(both,ref_allele > 0 & alt_allele > 0 & neighbor == j);
	temp4 <- subset(both,dm3_ref_ref_allele > 0 & dm3_alt_alt_allele > 0 & neighbor == j);
	
	medians[i,1] <- j;
	medians[i,2] <- median(log2(temp1$ref_allele.x/temp1$alt_allele.x));
	medians[i,3] <- median(log2(temp2$ref_allele.y/temp2$alt_allele.y));
	medians[i,4] <- median(log2(temp3$ref_allele/temp3$alt_allele));
	medians[i,5] <- median(log2(temp4$dm3_ref_ref_allele/temp4$dm3_alt_alt_allele));
	
	vars[i,1] <- j;
	vars[i,2] <- var(log2(temp1$ref_allele.x/temp1$alt_allele.x));
	vars[i,3] <- var(log2(temp2$ref_allele.y/temp2$alt_allele.y));
	vars[i,4] <- var(log2(temp3$ref_allele/temp3$alt_allele));
	vars[i,5] <- var(log2(temp4$dm3_ref_ref_allele/temp4$dm3_alt_alt_allele));
}

plot(medians[,1],medians[,2],type="o",col="black",xlim=c(0,13),ylim=c(-1,5),xaxt="n",yaxt="n",xlab="",ylab="");
par(new=TRUE);
plot(medians[,1],medians[,3],type="o",col="red",xlim=c(0,13),ylim=c(-1,5),xaxt="n",yaxt="n",xlab="",ylab="");
par(new=TRUE);
plot(medians[,1],medians[,4],type="o",col="green",xlim=c(0,13),ylim=c(-1,5),xlab="# neighboring SNPs",ylab="median log2(ASE)");

legend("topleft",legend=c("1 mismatch","2 mismatches","3 mismatches"),fill=c("black","red","green"));

# compare distributions of log2(ASE) for differing numbers of mismatches

total1 <- subset(complete_data,ref_allele.x > 0 & alt_allele.x > 0);
total2 <- subset(complete_data,ref_allele.y > 0 & alt_allele.y > 0);
total3 <- subset(complete_data,ref_allele > 0 & alt_allele > 0);

obj1 <- hist(log2(total1$ref_allele.x/total1$alt_allele.x),breaks=seq(-6,7,0.25));
obj2 <- hist(log2(total2$ref_allele.y/total2$alt_allele.y),breaks=seq(-6,7,0.25));
obj3 <- hist(log2(total3$ref_allele/total3$alt_allele),breaks=seq(-6,7,0.25));

mat <- cbind(obj1$counts/nrow(total1),obj2$counts/nrow(total2),obj3$counts/nrow(total3));

barplot(t(mat),beside=TRUE,names.arg=seq(-5.75,7,0.25),xlab="log2(ref/alt)",ylab="Proportion",main="SNP-based measurements of ASE",col=c("black","red","green"));

legend("topleft",legend=c("1 mismatch","2 mismatches","3 mismatches"),fill=c("black","red","green"));

#############################################

hist(log2(complete_data$ref_allele/complete_data$alt_allele));

complete_singleSNP <- subset(complete_data,complete_data$dist_left >= 50 & complete_data$dist_right >= 50);

hist(log2(complete_singleSNP$ref_allele/complete_singleSNP$alt_allele));

###############

data1 <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/bad_SNPs/ref_reads_bad_SNPs.dm3_ref.bowtie_v1_m1.SNP_ASE.txt",header=TRUE,sep="\t");

data2 <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/bad_SNPs/ref_reads_bad_SNPs.dm3_ref.bowtie_v1_m2.SNP_ASE.txt",header=TRUE,sep="\t");

data <- merge(data1,data2,by.x=c("chr","pos"),by.y=c("chr","pos"));

par(mfrow=c(1,2));
hist(log2(data$ref_allele.x/data$alt_allele.x),breaks=100,xlim=c(-6,6),ylim=c(0,150),xlab="log2(ref/alt)",main="1 mismatch, no multiple mapping");
hist(log2(data$ref_allele.y/data$alt_allele.y),breaks=100,xlim=c(-6,6),ylim=c(0,150),xlab="log2(ref/alt)",main="1 mismatch, maps at most 2 locations");

###################
# compare methods #
###################

##### compare single vs. multiple ASE measurements in exons

single <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/constExons_single_bp50_error0_tiled.dm3_ref.bowtie_v1_m1.exons.txt",header=TRUE,sep="\t");
multiple <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/constExons_single_bp50_error0_tiled.bowtie.exons.txt",header=TRUE,sep="\t");

single_multiple <- merge(single, multiple,by.x="gene_exon",by.y="gene_exon");

# > 0 ASE and AI

posASE_single <- subset(single_multiple,ref>0&alt>0&log2(ref/alt)!=0);
posASE_multiple <- subset(single_multiple,dm3_ref>0&dm3_alt>0&log2(dm3_ref/dm3_alt)!=0);

obj1 <- hist(log2(posASE_single$ref/posASE_single$alt),breaks=seq(-8,8,0.25));
obj2 <- hist(log2(posASE_multiple$dm3_ref/posASE_multiple$dm3_alt),breaks=seq(-8,8,0.25));
mat <- cbind(obj1$counts/nrow(posASE_single),obj2$counts/nrow(posASE_multiple));

barplot(t(mat),beside=TRUE,names.arg=seq(-7.75,8,0.25),xlab="log2(reference allele/alternative allele)",ylab="Proportion",main="Comparison of exon-based measurements\nof ASE between methods",col=c("black","gray"));

legend("topleft",legend=c("single (n = 31,334)","multiple (n = 674)"),fill=c("black","gray"));
legend("topright",legend="Exons with detectable ASE\nshowing imbalance");

pie(c(nrow(posASE_single),nrow(posASE_multiple)),labels=c("Single","Multiple"),col=c("black","gray"));

# > 0 ASE

posASE <- subset(single_multiple,ref>0&alt>0&dm3_ref>0&dm3_alt>0);

obj1 <- hist(log2(posASE$ref/posASE$alt),breaks=seq(-5,8,0.25));
obj2 <- hist(log2(posASE$dm3_ref/posASE$dm3_alt),breaks=seq(-5,8,0.25));
mat <- cbind(obj1$counts/nrow(posASE),obj2$counts/nrow(posASE));

barplot(t(mat),beside=TRUE,names.arg=seq(-4.75,8,0.25),xlab="",ylab="",main="",col=c("black","gray"));

legend("topleft",legend=c("single (n=42,809)","multiple (n=42,809)"),fill=c("black","gray"));
legend("topright",legend="Exons with detectable ASE");

##### compare single vs. multiple ASE measurements in SNPs

single <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/constExons_single_bp50_error0_tiled.dm3_ref.bowtie_v1_m1.SNP_ASE.txt",header=TRUE,sep="\t");
multiple <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/constExons_single_bp50_error0_tiled.bowtie.SNP_ASE.txt",header=TRUE,sep="\t");

single_multiple <- merge(single,multiple,by.x=c("chr","pos"),by.y=c("chr","pos"));

# > 0 ASE

posASE <- subset(single_multiple,ref_allele>0&alt_allele>0&dm3_ref_ref_allele>0&dm3_alt_alt_allele>0);

obj1 <- hist(log2(posASE$ref_allele/posASE$alt_allele),breaks=seq(-6,7,0.25));
obj2 <- hist(log2(posASE$dm3_ref_ref_allele/posASE$dm3_alt_alt_allele),breaks=seq(-6,7,0.25));
mat <- cbind(obj1$counts/nrow(posASE),obj2$counts/nrow(posASE));

barplot(t(mat),beside=TRUE,names.arg=seq(-5.75,7,0.25),xlab="log2(reference allele/alternative allele)",ylab="Proportion",main="SNP-based measurements of ASE",col=c("black","gray"));

legend("topleft",legend=c("single (n=243,718)","multiple (n=243,718)"),fill=c("black","gray"));
legend("topright",legend="SNPs with detectable ASE");

# > 0 ASE and AI

posASE_single <- subset(single_multiple,ref_allele>0&alt_allele>0&log2(ref_allele/alt_allele)!=0);
posASE_multiple <- subset(single_multiple,dm3_ref_ref_allele>0&dm3_alt_alt_allele>0&log2(dm3_ref_ref_allele/dm3_alt_alt_allele)!=0);

obj1 <- hist(log2(posASE_single$ref_allele/posASE_single$alt_allele),breaks=seq(-6,7,0.25));
obj2 <- hist(log2(posASE_multiple$dm3_ref_ref_allele/posASE_multiple$dm3_alt_alt_allele),breaks=seq(-6,7,0.25));
mat <- cbind(obj1$counts/nrow(posASE_single),obj2$counts/nrow(posASE_multiple));

barplot(t(mat),beside=TRUE,names.arg=seq(-5.75,7,0.25),xlab="log2(reference allele/alternative allele)",ylab="Proportion",main="Comparison of SNP-based measurements\nof ASE between methods",col=c("black","gray"));

legend("topleft",legend=c("single (n = 170,898)","multiple (n = 2,584)"),fill=c("black","gray"));
legend("topright",legend="SNPs with detectable ASE\nshowing imbalance");

pie(c(nrow(posASE_single),nrow(posASE_multiple)),labels=c("Single","Multiple"),col=c("black","gray"));

###############
# mappability #
###############

dm3_ref_exon_mappability <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/constExons.dm3_ref.l50_m0.mappability.txt",header=TRUE,sep="\t");

dm3_alt_exon_mappability <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/constExons.dm3_alt.l50_m0.mappability.txt",header=TRUE,sep="\t");

exon_mappability = merge(dm3_ref_exon_mappability,dm3_alt_exon_mappability,by.x="locus",by.y="locus");
dm3_ref_avg <- exon_mappability$sum.x/exon_mappability$length.x;
dm3_alt_avg <- exon_mappability$sum.y/exon_mappability$length.y;

mappability <- cbind(exon_mappability,dm3_ref_avg,dm3_alt_avg);

multiple <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/constExons_single_bp50_error0_tiled.bowtie.exons.txt",header=TRUE,sep="\t");

SNPsInExons <- read.table("/Users/kraigrs/Wittkopp/Simulations/SNPs_in_const.txt",header=FALSE,sep="\t");

temp1 <- merge(multiple,mappability,by.x="gene_exon",by.y="locus");
temp2 <- merge(temp1,SNPsInExons,by.x="gene_exon",by.y="V8");

weird <- subset(temp2,dm3_ref>0&dm3_alt>0&V7==0&log2(dm3_ref/dm3_alt)!=0)[,c(1:4,11,12,19)];

AI <- subset(temp2,dm3_ref>0&dm3_alt>0&V7>0&log2(dm3_ref/dm3_alt)!=0);
nrow(subset(AI,dm3_ref_avg != 1 | dm3_alt_avg != 1))/nrow(AI);
call <- rep("AI",nrow(AI));
AI <- cbind(AI,call);

ASE <- subset(temp2,dm3_ref>0&dm3_alt>0&V7>0&log2(dm3_ref/dm3_alt)==0);
nrow(subset(ASE,dm3_ref_avg != 1 | dm3_alt_avg != 1))/nrow(ASE);
call <- rep("ASE",nrow(ASE));
ASE <- cbind(ASE,call);

data <- rbind(AI,ASE);

boxplot(log2(dm3_ref_avg/dm3_alt_avg) ~ call, data=data,ylab="log2(mappability ratio)");

############

SNP_mappability <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/SNPs.dm3_ref.l50_m1.mappability.txt",header=TRUE,sep="\t");

single <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/constExons_single_bp50_error0_tiled.dm3_ref.bowtie_mm1.neighbor_SNPs.txt",sep="\t",header=TRUE);

data <- merge(single,SNP_mappability,by.x=c("chr","pos"),by.y=c("chr","position"));

AI <- subset(data,ref_allele > 0 & alt_allele > 0 & log2(ref_allele/alt_allele) != 0);
nrow(subset(AI,sum/length != 1))/nrow(AI);

ASE <- subset(data,ref_allele > 0 & alt_allele > 0 & log2(ref_allele/alt_allele) == 0);
nrow(subset(ASE,sum/length != 1))/nrow(ASE);

# SNPs with no neighbors and AI
temp <- subset(data,neighbor==0&log2(ref_allele/alt_allele)!=0);
nrow(subset(temp,sum/length != 1))/nrow(temp);

##################################
# plot mappability across genome #
##################################

library(lattice);

dm3_ref_exon_mappability <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/constExons.dm3_ref.l50_m0.mappability.txt",header=TRUE,sep="\t");

dm3_alt_exon_mappability <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/constExons.dm3_alt.l50_m0.mappability.txt",header=TRUE,sep="\t");

exon_mappability = merge(dm3_ref_exon_mappability,dm3_alt_exon_mappability,by.x="locus",by.y="locus");
dm3_ref_avg <- exon_mappability$sum.x/exon_mappability$length.x;
dm3_alt_avg <- exon_mappability$sum.y/exon_mappability$length.y;

mappability <- cbind(exon_mappability,dm3_ref_avg,dm3_alt_avg);

multiple <- read.table("/Users/kraigrs/Wittkopp/Simulations/tiled/constExons_single_bp50_error0_tiled.bowtie.exons.txt",header=TRUE,sep="\t");

SNPsInExons <- read.table("/Users/kraigrs/Wittkopp/Simulations/SNPs_in_const.txt",header=FALSE,sep="\t");

temp1 <- merge(multiple,mappability,by.x="gene_exon",by.y="locus");
temp2 <- merge(temp1,SNPsInExons,by.x="gene_exon",by.y="V8");

temp3 <- subset(temp2,dm3_ref > 0 & dm3_alt > 0);

xyplot(dm3_ref_avg/dm3_alt_avg ~ (V2+V3)/2 | V1, data = temp3,xlab="Midpoint of exon",ylab="ref_mappability/alt_mappability",main="Mappability across both allele-specific genomes",pch=19,cex=0.4,col=rgb(0,0,0,0.3));
xyplot(dm3_ref_avg/dm3_alt_avg ~ log2(dm3_ref/dm3_alt) | V1, data = temp3,xlab="log2(ref/alt) ASE",ylab="ref_mappability/alt_mappability",main="Concordance of mappability direction and bias in ASE",pch=19,cex=0.4,col=rgb(0,0,0,0.3));

# how many exons show opposite signs of mappability and bias? 0.2%, so 99.8% show same sign, awesome!
# also, of the 660 that do show differential mappability, 237 favor dm3_ref and 423 favor dm3_alt
sum( sign(log2(temp3$dm3_ref_avg/temp3$dm3_alt_avg)) != sign(log2(temp3$dm3_ref/temp3$dm3_alt)) );

# correlation? r^2 = 0.6
cor(log2(temp3$dm3_ref_avg/temp3$dm3_alt_avg),log2(temp3$dm3_ref/temp3$dm3_alt))^2;

