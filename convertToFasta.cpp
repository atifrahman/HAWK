#include <iostream>
#include <string>
#include <limits>
#include <vector>
#include <algorithm>
using namespace std;

#include <pthread.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#define MAX_REC_LEN 10240



int main(int argc, const char * argv[])
{
    
    	FILE *inFile=fopen("pvals_case_top_merged_sorted.txt","r");
    	FILE *outFile=fopen("case_kmers.fasta","w");
	FILE *totalFile=fopen("total_kmers.txt","r");
	

	char *line= new char[MAX_REC_LEN];
	char *line2= new char[MAX_REC_LEN];

	char *temp;

	double pval;

	long long int totalKmers;

	fscanf(totalFile,"%lld",&totalKmers);

    	int MAX_FILE_READ=MAX_REC_LEN/sizeof(line[0]);

 	while(fgets(line, MAX_FILE_READ, inFile)!=NULL)
       {
		temp=strtok(line,"\t\n ");
		pval=atof(temp);

		temp=strtok(NULL,"\t\n ");

		if(pval<0.05/totalKmers)
		{
			fprintf(outFile,">%e\n%s\n",pval,temp);
		}

	}



	fclose(inFile);
	fclose(outFile);
	fclose(totalFile);
	
 	inFile=fopen("pvals_control_top_merged_sorted.txt","r");
    	outFile=fopen("control_kmers.fasta","w");
	

 	while(fgets(line, MAX_FILE_READ, inFile)!=NULL)
       {
		temp=strtok(line,"\t\n ");
		pval=atof(temp);

		temp=strtok(NULL,"\t\n ");

		if(pval<0.05/totalKmers)
		{
			fprintf(outFile,">%e\n%s\n",pval,temp);
		}

	}



	fclose(inFile);
	fclose(outFile);
}