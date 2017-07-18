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

int main(int argc, const char * argv[])
{
 
	noCases=atoi(argv[1]);	
	noControls=atoi(argv[2]);    	

    char *temp;
    char kmerString[100];
    unsigned short count;
	
	int MAX_FILE_READ=10240;
	char *line=new char[MAX_FILE_READ+1];

	FILE * kmerFile=fopen("test.txt","r");
    	


	double SIG_LEVEL=0.05;


    FILE *kmerFileOut=fopen("kmerSummary.txt","w");
    
	long long int count1, count2;

	double meanCount1=0, meanCount2=0;

	double pVal;
	double meanpVal=0;

	int positive1=0, positive2=0;
	double meanPositive1=0, meanPositive2=0;


	long long int k=0;

	int val;

    	 while(fgets(line, MAX_FILE_READ, kmerFile)!=NULL)
        {
		k++;
            	temp=strtok(line,"\t\n ");
            	strcpy(kmerString,temp);
            
            	temp=strtok(NULL,"\t\n ");
//            	count1=atoll(temp);
		count1=0;		

		temp=strtok(NULL,"\t\n ");
//            	count2=atoll(temp);
		count2=0;  		
          			
		temp=strtok(NULL,"\t\n ");
            	pVal=atof(temp);

		positive1=0;
		positive2=0;

		for(int i=0;i<noCases;i++)
		{
			temp=strtok(NULL,"\t\n ");
	            	val=atoi(temp);
			count1+=val;
			if(val>0)
				positive1++;
		}

		for(int i=0;i<noControls;i++)
		{
			temp=strtok(NULL,"\t\n ");
	            	val=atoi(temp);
			count2+=val;
			if(val>0)
				positive2++;
		}
		meanCount1+=count1;
		meanCount2+=count2;
		meanPositive1+=positive1;
		meanPositive2+=positive2;
		meanpVal+=pVal;
	}


	meanCount1/=k;	
    	meanCount2/=k;
	meanPositive1/=k;
	meanPositive2/=k;
	meanpVal/=k;

	cout<<meanpVal<<",";
	cout<<meanCount1<<",";
	cout<<meanCount2<<",";
	cout<<meanPositive1<<",";
	cout<<meanPositive2<<endl;
	
	
    return 0;
}




