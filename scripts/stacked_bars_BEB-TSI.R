library(plotrix)
slices <- c(363224,96363,19584,28318,21798);
lgds <- c("SNPs", "Multimapped","Unmapped", "Indels", "Multiple SNPs/Structural" )
pct <- round(slices/sum(slices)*100,digits=1);
lbls <- paste("", pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
colors <- c("azure4","navy","gray","darkcyan","blue");

data=matrix(ncol=2,nrow=5);
texts=matrix(ncol=2,nrow=5);

data[,1]=pct;
texts[,1]=lbls;

slices <- c(281821,93223,18508,17905,50665); 
pct <- round(slices/sum(slices)*100,digits=1);
lbls <- paste("", pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="");data[,2]=pct;


texts[,2]=lbls;

coords=matrix(ncol=2,nrow=5);

for(i in 1:5)
{
	coords[i,2]=data[i,2]/2+sum(data[1:i-1,2]);
	coords[i,1]=data[i,1]/2+sum(data[1:i-1,1]);

}

bp<-barplot(data,width=c(1,1),col=colors,horiz=TRUE,space=c(5,5),names.arg=c("BEB","TSI"),las=1,xlim=c(0,100),ylim=c(0,20));

text(x=t(coords),y=bp-0.5,t(texts),cex=0.4,pos=1);

#pie3D(slices,labels=lbls,explode=0,main="BEB",col=colors);
legend("topright", lgds, cex=0.8, fill=colors)