#include <iostream>
#include <string>
#include <limits>
#include <vector>
#include <algorithm>
using namespace std;

//#include <pthread.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#define MAX_REC_LEN 10240
#define eps 1e-30


long long int findkcase(long long int totkmer){
	long long int m = totkmer;

	FILE *in=fopen("pvals_case_top_merged_sorted.txt","r");

	char *line= new char[MAX_REC_LEN];
	int MAX_FILE_READ=MAX_REC_LEN/sizeof(line[0]);
	
	long long int k = 1;
	long long int maxk = 0;
	while(fgets(line, MAX_FILE_READ, in)!=NULL){
       	char *temp=strtok(line,"\t\n ");
		double pval=atof(temp);
		
		double comp = 0.05*(((double)k)/(double)m);
		if(pval<comp || fabs(pval-comp)<eps){
			maxk = max(k,maxk);
		}
		k = k+1;
	}
	fclose(in);
	return maxk;
}

long long int findkcontrol(long long int totkmer){
	long long int m = totkmer;

	FILE *in=fopen("pvals_control_top_merged_sorted.txt","r");

	char *line= new char[MAX_REC_LEN];
	int MAX_FILE_READ=MAX_REC_LEN/sizeof(line[0]);
	
	long long int k = 1;
	long long int maxk = 0;
	while(fgets(line, MAX_FILE_READ, in)!=NULL){
       	char *temp=strtok(line,"\t\n ");
		double pval=atof(temp);
		
		double comp = 0.05*(((double)k)/(double)m);
		if(pval<comp || fabs(pval-comp)<eps){
			maxk = max(k,maxk);
		}
		k = k+1;
	}
	fclose(in);
	return maxk;
}

int main(int argc, const char * argv[])
{
    
    FILE *inFile=fopen("pvals_case_top_merged_sorted.txt","r");
    FILE *outFile=fopen("case_kmers_bh_correction.fasta","w");
	FILE *totalFile=fopen("total_kmers.txt","r");
	
	
	char *line= new char[MAX_REC_LEN];
	char *line2= new char[MAX_REC_LEN];

	char *temp;

	double pval;

	long long int totalKmers; 

	fscanf(totalFile,"%lld",&totalKmers);

	long long int kcase = findkcase(totalKmers);
	//cout<<kcase<<"\n";
    int MAX_FILE_READ=MAX_REC_LEN/sizeof(line[0]);

    long long int k = 1;
 	while(k <= kcase){
 		fgets(line, MAX_FILE_READ, inFile);
		temp=strtok(line,"\t\n ");
		pval=atof(temp);

		temp=strtok(NULL,"\t\n ");
		fprintf(outFile,">%e\n%s\n",pval,temp);
		
		k = k+1;
	}	



	fclose(inFile);
	fclose(outFile);
	//fclose(totalFile);
	
 	inFile=fopen("pvals_control_top_merged_sorted.txt","r");
    outFile=fopen("control_kmers_bh_correction.fasta","w");
	
    long long int kcontrol = findkcontrol(totalKmers);
    k = 1;
    //cout<<kcontrol<<"\n";
 	while(k <= kcontrol){
 		fgets(line, MAX_FILE_READ, inFile);
		temp=strtok(line,"\t\n ");
		pval=atof(temp);

		temp=strtok(NULL,"\t\n ");
		fprintf(outFile,">%e\n%s\n",pval,temp);
		
		k = k+1;

	}

	fclose(inFile);
	fclose(outFile);
	fclose(totalFile);
}
