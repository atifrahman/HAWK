p=0;
q=0;
coverage=5;
it=10000;

ps=c(1:20);
sigPvals_1=c(1:20);
sigPvals_2=c(1:20);

for(index in 1:20)
{

sig_1=0;
sig_2=0;

for(i in 1:it)
{
	x<-rmultinom(1,100,c(p^2,2*p*(1-p),(1-p)^2));

	case_0=x[1,1];
	case_1=x[2,1];
	case_2=x[3,1];

	p_case=(2*case_2+case_1)/(2*(case_2+case_1+case_0));

	y<-rmultinom(1,100,c(q^2,2*q*(1-q),(1-q)^2));


	control_0=y[1,1];
	control_1=y[2,1];
	control_2=y[3,1];
	
	p_control=(2*control_2+control_1)/(2*(control_2+control_1+control_0));
	
	total_0=case_0+control_0;
	total_1=case_1+control_1;
	total_2=case_2+control_2;
	
	p_total=(2*total_2+total_1)/(2*(total_2+total_1+total_0));
	
	l_alt=dmultinom(c(case_2,case_1,case_0),case_2+case_1+case_0,c(p_case^2,2*p_case*(1-p_case),(1-p_case)^2),log=TRUE)+dmultinom(c(control_2,control_1,control_0),control_2+control_1+control_0,c(p_control^2,2*p_control*(1-p_control),(1-p_control)^2),log=TRUE);
	
	l_null=dmultinom(c(case_2,case_1,case_0),case_2+case_1+case_0,c(p_total^2,2*p_total*(1-p_total),(1-p_total)^2),log=TRUE)+dmultinom(c(control_2,control_1,control_0),control_2+control_1+control_0,c(p_total^2,2*p_total*(1-p_total),(1-p_total)^2),log=TRUE);
	
	
	l_ratio=2*(l_alt-l_null);

	pval_1=(pchisq(l_ratio,1,lower.tail=FALSE));

	if(pval_1<0.05/39706715)
	{
		sig_1=sig_1+1;
	}

	control_count=0;

	if(y[2,1]>0)
	{
	for(j in 1:y[2,1])
	{
		val<-rpois(1,coverage/2);
		if(val>1)
		{
			control_count=control_count+val;
		}
	}
	}
	if(y[1,1]>0)
	{
	for(j in 1:y[1,1])
	{
		val<-rpois(1,coverage);
		if(val>1)
		{
			control_count=control_count+val;
		}
	}
	}

	case_count=0;
	if(x[2,1]>0)
	{
	for(j in 1:x[2,1])
	{
		val<-rpois(1,coverage/2);
		if(val>1)
		{
			case_count=case_count+val;
		}
	}
	}
	if(x[1,1]>0)
	{
	for(j in 1:x[1,1])
	{
		val<-rpois(1,coverage);
		if(val>1)
		{
			case_count=case_count+val;
		}
	}
	}
	mean_count=(case_count+control_count)/2;

	l_alt=dpois(case_count, case_count, log = TRUE)+dpois(control_count, control_count, log = TRUE);
	l_null=dpois(case_count, mean_count, log = TRUE)+dpois(control_count, mean_count, log = TRUE);
	
	l_ratio=2*(l_alt-l_null);

	pval_2=(pchisq(l_ratio,1,lower.tail=FALSE));

	if(pval_2<0.05/5598798965)
	{
		sig_2=sig_2+1;
	}
	
}

ps[index]=p;
sigPvals_1[index]=sig_1/it;
sigPvals_2[index]=sig_2/it;
p=p+0.01;

}

#plot(ps,sigPvals_1,type="l",col=1,xlab="Minor allele frequency",ylab="Significant fraction",xlim=c(0,0.2));
#legend(x=0.01,y=.95,legend="Multinomial",bty="n",lty=1,col=1);
#legend(x=0.01,y=.9,legend="Poisson (3x)",bty="n",lty=1,col=2);
#legend(x=0.01,y=.85,legend="Poisson (4x)",bty="n",lty=1,col=3);
#legend(x=0.01,y=.80,legend="Poisson (5x)",bty="n",lty=1,col=4);



lines(ps,sigPvals_2,col=4);

	