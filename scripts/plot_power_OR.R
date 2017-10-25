input=read.table("output_2_2.txt");

input[,2]=input[,2]/100;


plot(input[,1],input[,2],type="n",xlab="Number of samples",ylab="Power",ylim=c(0,1));

lw1 = glm(input[,2] ~ input[,1],family=quasibinomial,data=input)
lines(input[,1], lw1$fitted, col="green")

input=read.table("output_2_5.txt");
input[,2]=input[,2]/100;
lw1 = glm(input[,2] ~ input[,1],family=quasibinomial,data=input)
lines(input[,1], lw1$fitted, col="red")


input=read.table("output_2_10.txt");
input[,2]=input[,2]/100;
lw1 = glm(input[,2] ~ input[,1],family=quasibinomial,data=input)
lines(input[,1], lw1$fitted, col="blue")

legend(x=1500,y=0.5,legend="MAF=0.16",bty="n");
legend(x=1500,y=0.44,legend="OR=10",bty="n",col="blue",lty=1);
legend(x=1500,y=0.38,legend="OR=5",bty="n",col="red",lty=1);
legend(x=1500,y=0.32,legend="OR=2",bty="n",col="green",lty=1);