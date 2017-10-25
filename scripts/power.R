
totalKmer=10000000000;
alpha=0.5/totalKmer;

print(alpha);


k1=0;

powers<-rep(0,100);

for(d in 1:100)
{
	for(k2 in k1:100)
	{
		k_mean=(k1+k2)/2;

		lambda=dpois(k1,k1)*dpois(k2,k2)/(dpois(k1,k_mean)*dpois(k2,k_mean));

		prob=pchisq(2*log(lambda),1,lower.tail=FALSE);
	
		if(prob<=alpha)
		{
			powers[d]=(ppois(k2,d,lower.tail=FALSE));
			break;
		}

	}
}
plot(powers, type="l",col=8,xlab="Total k-mer coverage of cases",ylab="Power",xlim=c(0,60));
legend(x=40,y=0.5,legend="No of tests",bty="n");
legend(x=40,y=0.44,legend="10 billion",bty="n",col=8,lty=1);

for(i in 1:4)
{
alpha=alpha*10;
powers<-rep(0,100);

for(d in 1:100)
{
	for(k2 in k1:100)
	{
		k_mean=(k1+k2)/2;

		lambda=dpois(k1,k1)*dpois(k2,k2)/(dpois(k1,k_mean)*dpois(k2,k_mean));

		prob=pchisq(2*log(lambda),1,lower.tail=FALSE);
	
		if(prob<=alpha)
		{
			powers[d]=(ppois(k2,d,lower.tail=FALSE));
			break;
		}

	}
}

lines(powers,col=i);
}

legend(x=40,y=0.38,legend="1 billion",bty="n",col=1,lty=1);
legend(x=40,y=0.32,legend="100 million",bty="n",col=2,lty=1);
legend(x=40,y=0.26,legend="10 million",bty="n",col=3,lty=1);
legend(x=40,y=0.20,legend="1 million",bty="n",col=4,lty=1);



