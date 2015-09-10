mel_sec <- read.table("/Users/kraigrs/Desktop/mel_sec.txt",header=T,sep="\t");

coolonStevenson <- log2(mel_sec[,3]/mel_sec[,4]);

mcmanus <- log2(mel_sec[,6]/mel_sec[,7]);

r <- cor(coolonStevenson,mcmanus,method="pearson");

fit <- lm(coolonStevenson ~ mcmanus);
ci <- predict(fit,interval="confidence",type="response",level=0.95);

plot(mcmanus,coolonStevenson,
     pch=16,cex=0.5,col=rgb(0,0,0,0.4),
     xlim=c(-10,10),ylim=c(-10,10),
     main="log2(dm3/droSec1)");
    
for(i in 1:nrow(ci))   
{
  if(ci[i,1] < ci[i,2] || ci[i,1] > ci[i,3])
  {
    points(mcmanus[i],coolonStevenson[i],pch=16,col="red");
  }	
} 
    
abline(a=0,b=1);
abline(v=0,lty=2);
abline(h=0,lty=2);

model <- paste("lm: us = ",round(fit$coefficients[1],3)," + them*",round(fit$coefficients[2],3),sep="");
r2 = paste("r2 = ",round(r^2,3),sep="");

legend("topleft",legend=c(model,r2),bty="n");