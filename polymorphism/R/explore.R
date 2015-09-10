zhr_z30_5036 <- read.table("~/Wittkopp/polymorphism3/genesets/zhr_z30_5036_genes.txt",header=TRUE,sep="\t");
zhr_z30_7660 <- read.table("~/Wittkopp/polymorphism3/genesets/zhr_z30_7660_genes.txt",header=TRUE,sep="\t");

sim_sec_5036 <- read.table("~/Wittkopp/polymorphism3/genesets/sim_sec_5036_genes.txt",header=TRUE,sep="\t");
sim_sec_7660 <- read.table("~/Wittkopp/polymorphism3/genesets/sim_sec_7660_genes.txt",header=TRUE,sep="\t");

zhr_sim_5036 <- read.table("~/Wittkopp/polymorphism3/genesets/zhr_sim_5036_genes.txt",header=TRUE,sep="\t");
zhr_sim_7660 <- read.table("~/Wittkopp/polymorphism3/genesets/zhr_sim_7660_genes.txt",header=TRUE,sep="\t");

all_genes <- read.table("~/Wittkopp/polymorphism3/genesets/all_genes_chrs.txt",sep=" ");
colnames(all_genes) <- c("gene","mfbs","chrs");

#functions

gene_chr_props <- function(chrs,all)
{
	chrom_set <- data.frame(chrs);
	chrom_all <- data.frame(all);
	temp1 <- as.data.frame(table(chrom_set));
	temp2 <- as.data.frame(table(chrom_all));
	prop_set <- temp1$Freq/nrow(chrom_set);
	prop_all <- temp2$Freq/nrow(chrom_all);
	out <- cbind(temp1,temp2,prop_set,prop_all);
	return(out);
}

chr_seq_div <- function(merged)
{
	temp <- data.frame(levels(merged$chrs));
	colnames(temp) <- "chrs";
	
	diffs <- rep(0,nrow(temp));
	coding <- rep(0,nrow(temp));
	out <- cbind(temp,diffs,coding);

	for(i in 1:nrow(merged))
	{
		lev <- merged$chrs[i];
		j <- which(out$chrs == lev);
		out$diffs[j] <- out$diffs[j] + merged$diffs[i];
		out$coding[j] <- out$coding[j] + merged$coding[i];
	}
	seq_div <- out$diffs/out$coding;
	out <- cbind(out,seq_div);
	return(out);
}

#repeat this for each comparison to compare proportion of genes in set to all

data <- zhr_sim_7660;

chrs <- merge(data,all_genes,by="gene")$chrs;
all <- all_genes$chrs;
merged <- merge(data,all_genes,by="gene");

gene_chr_props(chrs,all);
chr_seq_div(merged);