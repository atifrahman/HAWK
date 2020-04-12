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

#define HASH_TABLE_LENGTH 50003

int kmerLength=31;

class Kmer
{
public:
	long long int val;
	double pVal_w_cor;
	double pVal_wo_cor;
	int countCase;
	int countControl;
	int positiveCase;
	int positiveControl;
};


class HashTable
{
public:
    vector <Kmer *> kmers[HASH_TABLE_LENGTH];
    HashTable();
    void insertKmer(long long int val, double pVal_w_cor, double pVal_wo_cor, int countCase, int countControl, int positiveCase, int positiveControl);
    Kmer * isKmerPresent(long long int val);	    
};

HashTable *ht;

long long int getInt(char *s)
{
	long long int val=0;
	int i=0;
	char ch;
	while(s[i])
	{
		val=val<<2;
		
		ch=s[i];
		if(ch=='A')
		{
			val=val|0;
		}
		else if(ch=='C')
		{
			val=val|1;
		}
		else if(ch=='G')
		{
			val=val|2;
		}
		else
		{
			val=val|3;
		}
		i++;
        
	}
	return val;
}

unsigned long int getHash(unsigned long long int key)
{
    /*
     key = (~key) + (key << 18); // key = (key << 18) - key - 1;
     key = key ^ (key >> 31);
     key = key * 21; // key = (key + (key << 2)) + (key << 4);
     key = key ^ (key >> 11);
     key = key + (key << 6);
     key = key ^ (key >> 22);
     */
	return (unsigned long int) key;
}

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


HashTable::HashTable()
{

}


void HashTable::insertKmer(long long int val,double pVal_w_cor, double pVal_wo_cor, int countCase, int countControl, int positiveCase, int positiveControl)
{
	unsigned long int index=getHash(val) % HASH_TABLE_LENGTH;
	int found=0;

    
	for(int i=0;i<kmers[index].size();i++)
	{
	 	if(val==kmers[index][i]->val)
        	{
			found=1;
			break;
		}
	}
	if(found==0)
	{

		Kmer* km=new Kmer;
		km->val=val;
        	km->pVal_wo_cor=pVal_wo_cor;
        	km->pVal_w_cor=pVal_w_cor;
        	km->countCase=countCase;
        	km->countControl=countControl;
        	km->positiveCase=positiveCase;
        	km->positiveControl=positiveControl;
		kmers[index].push_back(km);
       
	}
 
}

Kmer * HashTable::isKmerPresent(long long int val)
{
	unsigned long int index=getHash(val) % HASH_TABLE_LENGTH;
	int found=0;

    
	for(int i=0;i<kmers[index].size();i++)
	{
	 	if(val==kmers[index][i]->val)
        	{
			found=1;
			return kmers[index][i];
		}
	}
        
       return NULL;
 
}

int noCases, noControls;

vector<char*> contigs;
vector<char*> contigNames;
vector<unsigned long> contigLengths;
double *contigReadCounts;
long int noContigs;


int main(int argc, const char * argv[])
{

	noCases=atoi(argv[2]);	
	noControls=atoi(argv[3]); 	

	if(argc<2)
	{
		cout<<"Parameter missing: 1 for case stat, 0 for control stat"<<endl;
		exit(1);
	}
	int isCase=atoi(argv[1]);     

    	char *temp;
    	char kmerString[100];
	char reverseKmer[100];
    	unsigned short count;


	char fileName[1000];

	if(isCase==1)
	{
		strcpy(fileName,"pvals_case_top_merged_sorted.txt");
	}
	else
	{
		strcpy(fileName,"pvals_control_top_merged_sorted.txt");
	}
	
	int MAX_FILE_READ=10240;
	char *line=new char[MAX_FILE_READ+1];

	FILE * seedFile=fopen(argv[4],"r");
    	FILE * outFile=fopen("kmers.txt","w");

	FILE * statFile=fopen(fileName,"r");





		noContigs=0;
	
	long int contigLength=0;
	
	unsigned long read;
	char *contig;
	char *newcontig;
	char *contigName;
	long int bufferLength=1024;
    
	contigLength=0;
    
	long int tempContigLength=0;
    
	contig=new char[bufferLength];
	contig[0]='\0';
	
    // read contig file
    
	while(fgets(line, MAX_FILE_READ, seedFile)!=NULL)
	{
		if(line[0]==';')
		{
			continue;
		}
		else if(line[0]=='>')
		{
			contigName=new char[strlen(line)];
			strcpy(contigName,line+1);
			contigName[strlen(contigName)-1]='\0';
			contigNames.push_back(strtok(contigName," \t\n"));
			if(contigLength>0)
			{
				noContigs++;
				contig=(char *)realloc(contig, (contigLength+1)*sizeof(char));
				contigs.push_back(contig);
				contigLengths.push_back(contigLength);
				contigLength=0;
				bufferLength=1024;
				contig=new char[bufferLength];
				contig[0]='\0';
			}
		}
		else
		{
			read=strlen(line);
			tempContigLength=contigLength;
			if(read<MAX_FILE_READ-1)
			{
				contigLength+=(read-1);
			}
			else
			{
				contigLength+=MAX_FILE_READ-1;
				read++;
				
			}
			if(contigLength>bufferLength)
			{
				bufferLength=max(bufferLength*2,contigLength+1);
				newcontig=new char[bufferLength];
				strcpy(newcontig,contig);
				line[read-1]='\0';
				strcpy(newcontig+tempContigLength, line);
				delete []contig;
				contig=newcontig;
			}
			else
			{
				line[read-1]='\0';
				strcpy(contig+tempContigLength, line);
			}
            
		}
		
	}
    
    
    
	noContigs++;
	contig=(char *)realloc(contig, (contigLength+1)*sizeof(char));
	contigs.push_back(contig);
	contigLengths.push_back(contigLength);
	
    
 



	ht=new HashTable();

	long long int count1, count2;
	double pVal_wo_cor=0;
	double pVal_w_cor=0;


	int positive1=0, positive2=0;
	

	long long int k=0;

	int val;

    	 while(fgets(line, MAX_FILE_READ, statFile)!=NULL)
        {
		
            	temp=strtok(line,"\t\n ");
            	pVal_w_cor=atof(temp);

		temp=strtok(NULL,"\t\n ");
		strcpy(kmerString,temp);

		k++;

            
            	temp=strtok(NULL,"\t\n ");
		count1=0;		

		temp=strtok(NULL,"\t\n ");
		count2=0;  		
          			
		temp=strtok(NULL,"\t\n ");
            	pVal_wo_cor=atof(temp);

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

		ht->insertKmer(getInt(kmerString),pVal_w_cor,pVal_wo_cor,count1, count2, positive1, positive2);

	}


	Kmer *kmer;

	double meanCount1=0, meanCount2=0;

	double meanpVal_wo_cor=0;
	double meanpVal_w_cor=0;
	double meanPositive1=0, meanPositive2=0;


	 for(int i=0;i<noContigs;i++)	
//    	 while(fgets(line, MAX_FILE_READ, seedFile)!=NULL)
        {
		
		line=contigs[i];

		meanpVal_wo_cor=0;
		meanpVal_w_cor=0;
		meanCount1=0;
		meanCount2=0;
		meanPositive1=0;
		meanPositive2=0;
		k=0;
	
		for(int i=0;i<strlen(line)-kmerLength;i++)
		{
			strncpy(kmerString, line+i, kmerLength);
			kmer=ht->isKmerPresent(getInt(kmerString));

			if(kmer==NULL)
				continue;

			k++;
			meanpVal_wo_cor+=kmer->pVal_wo_cor;
			meanpVal_w_cor+=kmer->pVal_w_cor;;
			meanCount1+=kmer->countCase;
			meanCount2+=kmer->countControl;
			meanPositive1+=kmer->positiveCase;
			meanPositive2+=kmer->positiveControl;
			

			
			reverse(reverseKmer,kmerString);

			kmer=ht->isKmerPresent(getInt(reverseKmer));
			if(kmer==NULL)
				continue;	

			k++;
			meanpVal_wo_cor+=kmer->pVal_wo_cor;
			meanpVal_w_cor+=kmer->pVal_w_cor;;
			meanCount1+=kmer->countCase;
			meanCount2+=kmer->countControl;
			meanPositive1+=kmer->positiveCase;
			meanPositive2+=kmer->positiveControl;

	
       	}     
				
		if(k>0)
		{			
			meanCount1/=k;	
    			meanCount2/=k;
			meanPositive1/=k;
			meanPositive2/=k;
			meanpVal_w_cor/=k;
			meanpVal_wo_cor/=k;

			cout<<meanpVal_w_cor<<",";
			cout<<meanpVal_wo_cor<<",";
			cout<<meanCount1<<",";
			cout<<meanCount2<<",";
			cout<<meanPositive1<<",";
			cout<<meanPositive2<<endl;

		}
		else
		{
			cout<<-1<<",";
			cout<<-1<<",";
			cout<<-1<<",";
			cout<<-1<<",";
			cout<<-1<<",";
			cout<<-1<<endl;

				
		}

	}

	return 0;
}
