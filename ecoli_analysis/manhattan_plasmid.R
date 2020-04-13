sam=read.table("case_kmers_r_plasmid.sam",fill=TRUE);

length=107229;

plot(sam[,4],-log10(sam[,1]),xlim=c(0,length),ylim=c(0,30),pch=19,col="navy",xaxt="n",xlab="Plasmid pKBN10P04869A",ylab=expression(paste(-log[10], "(p-value)")));



axis(1,at=c(0,25000,50000,75000,100000), labels=c(0,25000,50000,75000,"100000"));

abline(v=14596,col="azure4",lwd=1);

abline(v=38477,col="azure4",lwd=1);

text("blaTEM",x=14596,y=28,pos=4);

text("blaTEM",x=38477,y=28,pos=4);