BEB<-c(529287,413340+96363,7560+7894,15337+14630,2599+2862);
TSI<-c(462122,350391+93223,12542+10972,19319+19113,3483+5194);

BEBtext=BEB/529287*100;
TSItext=TSI/462122*100;


counts<-matrix(c(BEB,TSI),nrow=5,ncol=2);
texts<-matrix(c(BEBtext,TSItext),nrow=5,ncol=2);

colors <- c("navy","azure4");

bp<-barplot(t(counts),beside=TRUE,legend=c("BEB","TSI"),names.arg=c("Total","hg38","RefSeq","Exons","Coding"),col=colors,ylab="Number of sequences",ylim=c(0,550000),yaxt='n');
text(bp,t(counts),paste(" ",round(t(texts),1),"%"),pos=3,cex=0.6);

axis(side=2, at=c(0,100000,200000,300000,400000,500000), labels=c("0","100000","200000","300000","400000","500000"));
