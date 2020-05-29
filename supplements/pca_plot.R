#library(scatterplot3d)


PC1=1;
PC2=2;


# read in the vcf2eigen.pca.evec file
df <- read.table("gwas_eigenstrat.evec");

df[,c(12)] <- sapply(df[,c(12)],as.character) 

for(i in 1:nrow(df))
{
	if(df[i,12]=="Case")
		df[i,12]=1;
	if(df[i,12]=="Control")
		df[i,12]=2;
}

names(df) <- c("SampleID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "POPULATION")

setEPS()
postscript('pca_plot.eps')
plot(df[,PC1+1], df[,PC2+1],  col = df[,12], pch=19, main="PCA plot",xlab="PC1",ylab="PC2");

legend("topright", pch=19, col=c(1,2), c("Case", "Control"));
dev.off()



