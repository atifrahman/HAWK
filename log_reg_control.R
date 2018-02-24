library(foreach)
library(doParallel)

cl<-makeCluster(1)
registerDoParallel(cl)


con = file("control_out_w_bonf_top.kmerDiff", "r")

Z <- read.table("pcs.evec",header=FALSE);

input <- read.table("gwas_eigenstratX.ind");

Y <- matrix(nrow=nrow(Z),ncol=1);
cov <- matrix(nrow=nrow(Z),ncol=1);

totals <- read.table("total_kmer_counts.txt");

counts<- vector(length=length(Y));


for(i in 1:length(Y))
{
	if(input[i,3]=="Case")
	{
		Y[i,1]=1;
	}
	else if(input[i,3]=="Control")
	{
		Y[i,1]=0;
	}
	else
	{
		Y[i,1]=input[i,3];
	}
}

print(Y);

if (FALSE) { 
input <- read.table("covariates.txt");
for(i in 1:length(cov))
{
	if(input[i,1]=="M")
	{
		cov[i,1]=0;
	}
	else
	{
		cov[i,1]=1;
	}
}
}

totalKmers <- read.table("total_kmers.txt",header=FALSE);

n<-nrow(Z);

cat("", file = "pvals_control_top.txt")

CHUNK_SIZE=10000;

ptm <- proc.time()


while(TRUE)
{
kmercounts <- read.table(con,nrow=CHUNK_SIZE);

nr=nrow(kmercounts);


ls<-foreach(j=icount(nr), .combine=cbind) %dopar%
{


	for(i in 1:length(Y))
	{
		counts[i]=kmercounts[j,4+i]/totals[i,1];
	}

#	model1<-glm(formula = Y ~ counts, family = binomial(link = "logit"));
	

	model2<-glm(formula = Y ~ Z[,1]+Z[,2]+totals[,1]+counts, family = binomial(link = "logit"));

	#summary(model1);

#	v1<-anova(model1, test="Chisq");

	#summary(model2);

	v2<-anova(model2, test="Chisq");

#	rbind(v1$'Pr(>Chi)'[2],v2$'Pr(>Chi)'[13]);

	v2$'Pr(>Chi)'[5];
}

write(ls,file='pvals_control_top.txt',ncolumns=1,append=TRUE,sep='\t');

if(nr<CHUNK_SIZE)
{
	break;
}


	
}


close(con); 

print(proc.time() - ptm)


stopCluster(cl)

