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

char * toLower(char *s)
{
	int i=0;
	while(s[i]!=0)
	{
		s[i]=tolower(s[i]);
		i++;
	}
	return s;

}

int main(int argc, const char * argv[])
{
    
    
 
    	FILE *kmerFile=fopen("sorted_files.txt","r");
    	FILE *totalFile=fopen("total_kmer_counts.txt","r");
	FILE *infoFile=fopen("gwas_info.txt","r");

	FILE *caseKmerFile=fopen("case_sorted_files.txt","w");
    	FILE *caseTotalFile=fopen("case_total_kmers.txt","w");
	FILE *caseIndFile=fopen("case.ind","w");


	FILE *controlKmerFile=fopen("control_sorted_files.txt","w");
    	FILE *controlTotalFile=fopen("control_total_kmers.txt","w");
	FILE *controlIndFile=fopen("control.ind","w");

	

	char *line= new char[MAX_REC_LEN];
	char *line2= new char[MAX_REC_LEN];
	char *line3= new char[MAX_REC_LEN];


    	int MAX_FILE_READ=MAX_REC_LEN/sizeof(line[0]);

    	char *temp;

	int count,flag=0,n;

	while(fgets(line, MAX_FILE_READ, infoFile)!=NULL)
       {

		strcpy(line2,line);

		temp=strtok(line,"\t\n ");
		temp=strtok(NULL,"\t\n ");
		temp=strtok(NULL,"\t\n ");

		if(strncmp(toLower(temp),"case",4)==0 || strncmp(toLower(temp),"1",1)==0)
		{
			fputs(line2,caseIndFile);


			fgets(line2, MAX_FILE_READ, kmerFile);
			fputs(line2,caseKmerFile);

			fgets(line2, MAX_FILE_READ, totalFile);
			fputs(line2,caseTotalFile);

		}
		else if(strncmp(toLower(temp),"control",7)==0  || strncmp(toLower(temp),"0",1)==0)
		{
			fputs(line2,controlIndFile);


			fgets(line2, MAX_FILE_READ, kmerFile);
			fputs(line2,controlKmerFile);

			fgets(line2, MAX_FILE_READ, totalFile);
			fputs(line2,controlTotalFile);



		}


	}



	fclose(kmerFile);
	fclose(totalFile);
	fclose(infoFile);

	
	fclose(caseKmerFile);
	fclose(caseTotalFile);
	fclose(caseIndFile);
	fclose(controlKmerFile);
	fclose(controlTotalFile);
	fclose(controlIndFile);
	


}