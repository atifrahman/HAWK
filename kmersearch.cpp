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

int kmerLength=31;

void reverse(char *reverse, char *read)
{
    
	char ch='A';
	int readLength=strlen(read);
	for(int i=1; i<=readLength;i++)
	{
		ch=read[i-1];
		if(ch=='A')
			reverse[readLength-i]='T';
		else if(ch=='C')
			reverse[readLength-i]='G';
		else if(ch=='G')
			reverse[readLength-i]='C';
		else if(ch=='T')
			reverse[readLength-i]='A';
		else
			reverse[readLength-i]='N';
	}
	reverse[readLength]='\0';
    
}


int main(int argc, const char * argv[])
{
     
    char *temp;
    char kmerString[100];
	char reverseKmer[100];
    unsigned short count;
	
	int MAX_FILE_READ=10240;
	char *line=new char[MAX_FILE_READ+1];

	FILE * seedFile=fopen("seed.txt","r");
    	FILE * outFile=fopen("kmers.txt","w");



    	 while(fgets(line, MAX_FILE_READ, seedFile)!=NULL)
        {

		for(int i=0;i<strlen(line)-kmerLength+1;i++)
		{
			strncpy(kmerString, line+i, kmerLength);
			fprintf(outFile,"%s\n",kmerString);
			reverse(reverseKmer,kmerString);
			fprintf(outFile,"%s\n",reverseKmer);
			
	
       	}     
		
	}
    
    return 0;
}




