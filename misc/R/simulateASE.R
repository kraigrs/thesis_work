n <- 10000;    # number of data points (could be number of genes)
lambda <- 100; # average coverage; parameter for Poisson draws

binoms <- NULL;
for(i in 1:n)
{
	tot <- rpois(1,lambda); # draw from Poisson(lambda) to get coverage
	binoms <- rbind(binoms,c(rbinom(1,tot,0.5),tot)); # simulate successes from tot with p=0.5
}

pvals <- NULL;
for(i in 1:nrow(binoms))
{
	pvals <- c(pvals,binom.test(binoms[i,1],binoms[i,2],p=0.5,alternative="two.sided",conf.level=0.95)$p.value);
}

qvals <- p.adjust(pvals,method="fdr");
x <- runif(n); # to be used in the generation of qq-plots

# qq-plots

plot(sort(x),sort(pvals),xlim=c(0,1),ylim=c(0,1),xlab="Theoretical quantiles",ylab="p-values from binom.test for ASE");

# histogram

hist(pvals,breaks=50);