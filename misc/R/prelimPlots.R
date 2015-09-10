# Make plots for prelim paper

mel <- rnorm(1000,mean=0,sd=1);
yak <- rnorm(1000,mean=-1.5,sd=1);

melvar = mean(mel);
yakvar = mean(yak)

xlimits <- c(min(mel,yak),max(mel,yak));
ylimits <- c(0,max(max(hist(mel,plot=F)$density,hist(yak,plot=F)$density)));

brk <- seq(floor(min(mel,yak)),ceiling(max(mel,yak)),0.5);

melColor <- rgb(t(col2rgb("blue")),alpha=50,maxColorValue=255);
yakColor <- rgb(t(col2rgb("yellow")),alpha=50,maxColorValue=255);

plot(density(mel),xlim=xlimits,ylim=ylimits,xlab="",yaxt="n",ylab="",main="");
par(new=TRUE);
hist(mel,freq=FALSE,breaks=brk,col=melColor,xlim=xlimits,ylim=ylimits,xlab="",yaxt="n",ylab="",main="");
par(new=TRUE);
plot(density(yak),xlim=xlimits,ylim=ylimits,xlab="",yaxt="n",ylab="",main="");
par(new=TRUE);
hist(yak,freq=FALSE,breaks=brk,col=yakColor,xlim=xlimits,ylim=ylimits,xlab="log2(s1/s2)",yaxt="n",ylab="Frequency",main="");
legend("topright",legend=c("corrected","uncorrected"),fill=c(melColor,yakColor));


