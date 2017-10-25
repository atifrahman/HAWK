YRI<-c(2970929,2523855+377823,40419+36035,82469+65190,7711+14803);
TSI<-c(1865285,1540647+280276,42190+36226,68094+64091,10087+19421);

YRItext=YRI/2970929*100;
TSItext=TSI/1865285*100;


counts<-matrix(c(YRI,TSI),nrow=5,ncol=2);
texts<-matrix(c(YRItext,TSItext),nrow=5,ncol=2);

colors <- c("navy","azure4");

bp<-barplot(t(counts),beside=TRUE,legend=c("YRI","TSI"),names.arg=c("Total","hg38","RefSeq","Exons","Coding"),col=colors,ylab="Number of sequences",ylim=c(0,3200000));
text(bp,t(counts),paste(" ",round(t(texts),1),"%"),pos=3,cex=0.6);