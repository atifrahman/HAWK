noSamples=150;

for(noSamples in 100:101)
{
S_R=1/3;

noCases=noSamples*S_R;
noControls=noSamples-noCases;

MAF=0.16;

noExposed=noSamples*MAF;
noNotExposed=noSamples-noExposed;

oddsRatio=5;

a=oddsRatio-1;
b=oddsRatio*(-noCases-noExposed)-(noSamples-noCases-noExposed);
c=noCases*noExposed*oddsRatio;

x1=(-b+sqrt(b^2-4*a*c))/(2*a);
x2=(-b-sqrt(b^2-4*a*c))/(2*a);

D_E=0;

if ((x1<0 && x2<0) || (x1>noCases && x2>noCases))
{
	print("error");
}else if (x2<0)
{
	D_E=x1;
}else if (x1>noCases)
{
	D_E=x2;
}else
{
	D_E=x1;
}

P_e_case=D_E/noCases;

P_e_control=(noExposed-D_E)/noControls;

noTrials=100;
noRejected=0;

for(it in 1:noTrials)
{
cases=rbinom(noCases,1,P_e_case);
controls=rbinom(noControls,1,P_e_control);

cov=5;

caseCounts=rpois(noCases,cov);
controlCounts=rpois(noControls,cov);

for(i in 1:noCases)
{
	caseCounts[i]=caseCounts[i]*cases[i];
}

for(i in 1:noControls)
{
	controlCounts[i]=controlCounts[i]*controls[i];
}

k_case=sum(caseCounts)/noCases;
k_control=sum(controlCounts)/noControls;

k_all=(sum(caseCounts)+sum(controlCounts))/noSamples;

likelihood_alt=0;
likelihood_null=0;

for(i in 1:noCases)
{
	likelihood_alt=likelihood_alt+ppois(caseCounts[i],k_case,log.p = TRUE);
	likelihood_null=likelihood_null+ppois(caseCounts[i],k_all,log.p = TRUE);
}
for(i in 1:noControls)
{
	likelihood_alt=likelihood_alt+ppois(controlCounts[i],k_control,log.p = TRUE);
	likelihood_null=likelihood_null+ppois(controlCounts[i],k_all,log.p = TRUE);
}


lambda=2*(likelihood_alt-likelihood_null);

prob=pchisq(lambda,1,lower.tail=FALSE);

if(prob<0.05/10000)
{
	noRejected=noRejected+1;
}
}
print(noRejected);
cat(noSamples," ",noRejected,"\n", file="output.txt",append=TRUE);
}