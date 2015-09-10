##### compare significance of binomial exact tests and G-tests #####

# define parameters

N <- 100; # number of expected reads to overlap SNPs (this is really only an assumption for DNA)
p <- 0.5; # expected
SNPs <- 50000; # number of SNPs being tested

df <- N - 1;

# generate numbers of reference reads with equal probability at SNPs

ref_reads <- rbinom(SNPs,N,p);

# perform tests and compare results

binom_tests_pvals <- mat.or.vec(SNPs,1);
G_tests_pvals <- mat.or.vec(SNPs,1);

for(i in 1:SNPs)
{
	temp1 <- binom.test(ref_reads[i],N,p=p,alternative="two.sided",conf.level=0.95);
	G <- 2*(ref_reads[i]*log(ref_reads[i]/(N/2))+(N-ref_reads[i])*log((N-ref_reads[i])/(N/2)));
	temp2 <- pchisq(G,df=df,lower.tail=FALSE);
	
	binom_tests_pvals[i] <- temp1$p.value;
	G_tests_pvals[i] <- temp2;
}

binom_tests_qvals <- p.adjust(binom_tests_pvals,method="fdr");
G_tests_qvals <- p.adjust(G_tests_pvals,method="fdr");

# make QQ plot comparing these tests

binom_sorted <- sort(binom_tests_pvals);
G_sorted <- sort(G_tests_pvals);

plot(-log10(binom_sorted),-log10(G_sorted),
	xlab="-log10(p-value from binomial exact test) quantiles",
	ylab="-log10(p-value from G-test) quantiles",
	main="Compare tests of differential ASE",
	xlim=c(0,0.4),
	ylim=c(0,0.4)
	);
	
abline(a=0,b=1,col="red",lty=2);
	
	
	
	