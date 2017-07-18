//
//  main.cpp
//  kmer
//

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

int noCases, noControls;

int kmerLength=31;

int NUM_THREADS=16;


int main(int argc, const char * argv[])
{
 

    
    char *temp;
    char kmerString[100];
    unsigned short count;
	
	int MAX_FILE_READ=10240;
	char *line=new char[MAX_FILE_READ+1];

	long long int TOTAL_KMER;
	double SIG_LEVEL=0.05;

	FILE * kmerFile=fopen("total_kmers.txt","r");
	
	fscanf(kmerFile,"%lld",&TOTAL_KMER);

	FILE * caseFile=fopen("case_out_wo_bonf.kmerDiff","r");
    	FILE * controlFile=fopen("control_out_wo_bonf.kmerDiff","r");


    FILE *caseFileOut=fopen("case_out_sig.fasta","w");
    FILE *controlFileOut=fopen("control_out_sig.fasta","w");

	int count1, count2;
	double pVal;

	long long int i=0;

    	 while(fgets(line, MAX_FILE_READ, caseFile)!=NULL)
        {
		temp=strtok(line,"\t\n ");
            	strcpy(kmerString,temp);
            
            	temp=strtok(NULL,"\t\n ");
            	count1=atoi(temp);

		temp=strtok(NULL,"\t\n ");
            	count2=atoi(temp);
            			
		temp=strtok(NULL,"\t\n ");
            	pVal=atof(temp);


		if(pVal<SIG_LEVEL/TOTAL_KMER)
		{
			fprintf(caseFileOut,">%lld\n",i);
			fprintf(caseFileOut,"%s\n",kmerString);
               
			i++;
		}
	}

	i=0;

	while(fgets(line, MAX_FILE_READ, controlFile)!=NULL)
        {
            	temp=strtok(line,"\t\n ");
            	strcpy(kmerString,temp);
            
            	temp=strtok(NULL,"\t\n ");
            	count1=atoi(temp);

		temp=strtok(NULL,"\t\n ");
            	count2=atoi(temp);
            			
		temp=strtok(NULL,"\t\n ");
            	pVal=atof(temp);

		if(pVal<SIG_LEVEL/TOTAL_KMER)
		{
			fprintf(controlFileOut,">%lld\n",i);
		
			fprintf(controlFileOut,"%s\n",kmerString);
               
			i++;
		}
	}

    
    return 0;
}




