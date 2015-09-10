####################################################
# script to determine imprinting and misexpression #
####################################################

FET_misexpression <- function(locus)
{
	table <- matrix(c(locus[1],locus[3]-locus[1],locus[2],locus[4]-locus[2]),
					nrow=2,
					dimnames=list(c("gene","not gene"),c("sample1","sample2")));
	test <- fisher.test(table,or=1,alternative="two.sided",conf.level=0.95);
	result <- c(test$estimate,test$conf.int[1],test$conf.int[2],test$p.value);
	return(result);
}

FET_imprinting <- function(locus)
{
	table <- matrix(c(locus[1],locus[3],locus[2],locus[4]),
					nrow=2,
					dimnames=list(c("H1","H2"),c("Allele1","Allele2")));
	test <- fisher.test(table,or=1,alternative="two.sided",conf.level=0.95);
	result <- c(test$estimate,test$conf.int[1],test$conf.int[2],test$p.value);
	return(result);
}

boot <- function(mat)
{
	set1 <- sample(seq(1,nrow(mat),1),nrow(mat),replace=TRUE);
	set2 <- sample(seq(1,nrow(mat),1),nrow(mat),replace=TRUE);
	val1 <- 1-cor(log10(mat$Dpse_TL.x[set1]),log10(mat$Dbog_Toro1.x[set1]),use="pairwise.complete.obs",method="spearman");
	val2 <- 1-cor(log10(mat$Dpse_TL.y[set2]),log10(mat$Dbog_Toro1.y[set2]),use="pairwise.complete.obs",method="spearman");
	return(c(val1,val2));
}

################
# read in data #
################

reg_div1 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_ovaries_reg_div.txt",header=TRUE,sep="\t");
MOI1 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H6_ovaries_MOI.txt",header=TRUE,sep="\t");
data1 <- merge(reg_div1,MOI1,by.x=1,by.y=1);

reg_div2 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H5_ovaries_reg_div.txt",header=TRUE,sep="\t");
MOI2 <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/data/TL_Toro1_H5_ovaries_MOI.txt",header=TRUE,sep="\t");
data2 <- merge(reg_div2,MOI2,by.x=1,by.y=1);

data <- merge(data1,data2,by.x="gene",by.y="gene");

annotation <- read.table("/Users/kraigrs/Wittkopp/Machado_Dpse/genomes/dpse_r3.1.genes.bed",header=FALSE,sep="\t");
colnames(annotation) <- c("chromosome","start","stop","gene","empty","strand");

data_annotation <- merge(data,annotation,by.x="gene",by.y="gene");

##################
# female samples #
##################

# imprinting
imprinting_results <- t(apply(as.matrix(data_annotation[,c(4,5,29,30)]),1,FET_imprinting));
imprinting_qval <- p.adjust(imprinting_results[,4],method="fdr");
data_imprinting <- cbind(data_annotation,imprinting_qval);
sum(data_imprinting$imprinting_qval<0.05);

# hybrid misexpression
misexpression_results <- t(apply(as.matrix(cbind(data_annotation[,c(25,50)],rep(sum(data_annotation[,25],nrow(data_annotation))),rep(sum(data_annotation[,50],nrow(data_annotation))))),1,FET_misexpression));
misexpression_qval <- p.adjust(misexpression_results[,4],method="fdr");
data_misexpression <- cbind(data_annotation,misexpression_qval);
sum(data_misexpression$misexpression_qval<0.05);

################
# male samples #
################

data_no_X <- data_annotation[grep("^[^X]",data_annotation$chromosome,perl=TRUE),];

# imprinting
imprinting_results <- t(apply(as.matrix(data_no_X[,c(4,5,29,30)]),1,FET_imprinting));
imprinting_qval <- p.adjust(imprinting_results[,4],method="fdr");
data_imprinting <- cbind(data_no_X,imprinting_qval);
sum(data_imprinting$imprinting_qval<0.05);

# hybrid misexpression
misexpression_results <- t(apply(as.matrix(cbind(data_no_X[,c(25,50)],rep(sum(data_no_X[,25],nrow(data_no_X))),rep(sum(data_no_X[,50],nrow(data_no_X))))),1,FET_misexpression));
misexpression_qval <- p.adjust(misexpression_results[,4],method="fdr");
data_misexpression <- cbind(data_no_X,misexpression_qval);
sum(data_misexpression$misexpression_qval<0.05);




















