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
#include "kmer.h"

#include "specialfunctions.h"

#define VERSION "0.9.8-beta"

#define MAX_REC_LEN 10240

int noCases, noControls;

int kmerLength=31;

int NUM_THREADS=16;

pthread_mutex_t readFile_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t caseOutFile_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t controlOutFile_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t hashTable_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t eigenFile_mutex = PTHREAD_MUTEX_INITIALIZER;


FILE *eigenGenoFile;
FILE *eigenSNPFile;



#pragma pack(push, 1)
class Kmer
{
public:
	long long int kmer;
	unsigned short * caseCounts;
    	unsigned short * controlCounts;
	double pVal;
	double meanCase;
	double meanControl;
    char significanceType;
	char forPCA;
	char forPval;


    
    Kmer(int noCases, int noControls);
    ~Kmer();
    void show();
	void freeMemory();
};
#pragma pack(pop)

class KeyVal
{
public:
	long long int val;
	int count;
};

class HashTable
{
public:
    long long int totalKmers;
    long long int totalTests;
    vector <Kmer *> kmers[HASH_TABLE_LENGTH];
	pthread_mutex_t hashTableBuckets_mutex[HASH_TABLE_LENGTH];

	long long int totalKmers_Bucket[HASH_TABLE_LENGTH];	

	
	long long int * totalKmerCountsCase;	
	long long int * totalKmerCountsControl;
    HashTable();
    void insertKmer(long long int val, int count, int isCase, int sampleNo);
    void show();
    void computePvalues();
    void computeLikelihoodRatios();	
    void dumpKmers(double sigLevel);
    
};

HashTable *ht;

class Factorials
{
public:
    double *factorials;
    Factorials(int max);
    double getFactorial(int number);
};

Factorials *f;

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

char * getKmer(long long int val, char *kmer, int kmerLength)
{
    
	int temp=0;
	for(int i=kmerLength-1;i>=0;i--)
	{
		temp=val&3;
		val=val>>2;
		if(temp==0)
			kmer[i]='A';
		else if(temp==1)
			kmer[i]='C';
		else if(temp==2)
			kmer[i]='G';
		else if(temp==3)
			kmer[i]='T';
	}
	kmer[kmerLength]='\0';
	return kmer;
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

double getPvalue(int a, int b, int c, int d)
{
    int a_b=a+b;
    int a_c=a+c;
    int b_d=b+d;
    int c_d=c+d;
    int n=a_b+c_d;
    
    double pVal=0;
    double p=a_b/(double)n;
    double temp;
    
    if(a/(double)a_c<p)
    {
        for(int i=a;i>=0;i--)
        {
            temp=f->getFactorial(i)*f->getFactorial(a_b-i)*f->getFactorial(a_c-i)*f->getFactorial(b_d-a_b+i);
            pVal+=1/temp;
            
        }
    }
    else
    {
        for(int i=a;i<=a_b;i++)
        {
            temp=f->getFactorial(i)*f->getFactorial(a_b-i)*f->getFactorial(a_c-i)*f->getFactorial(b_d-a_b+i);
            pVal+=1/temp;
            
        }
        
    }
    pVal/=f->getFactorial(n);
    temp=(f->getFactorial(a_b)*f->getFactorial(c_d)*f->getFactorial(a_c)*f->getFactorial(b_d));
    pVal*=temp;
    
    
    return pVal;
}

double getPvalue_chi(int a, int b, int c, int d)
{
    double n=a+b+c+d;
    double a_E=(a+b)*(a+c)/n;
    double b_E=(a+b)*(b+d)/n;;
    double c_E=(a+c)*(c+d)/n;
    double d_E=(b+d)*(c+d)/n;
    
    double statistic=(a-a_E)*(a-a_E)/a_E;
    statistic+=(b-b_E)*(b-b_E)/b_E;
    statistic+=(c-c_E)*(c-c_E)/c_E;
    statistic+=(d-d_E)*(d-d_E)/d_E;
    
    double pVal=alglib::chisquarecdistribution(1, statistic);
    
    
    return pVal;
}

double logFactorial( int k )
{
   double result = 0;
   for(; k > 1; k--)
   {
       result += log(k);
   }

   return result;
}


double poissonProb(int k, double lambda)
{  
	if(lambda<=0)
		return 0;
	if(k<0)
		k=0;
	
//	return	exp( -lambda + (k*log(lambda)) - logFactorial(k) );
  	return	( -lambda + (k*log(lambda)) - logFactorial(k) );
  

//    return alglib::poissondistribution(count, rate)-alglib::poissondistribution(count-1, rate);
    
}


Kmer::Kmer (int noCases, int noControls)
{
	caseCounts=new unsigned short[noCases];
    for(int i=0;i<noCases;i++)
        caseCounts[i]=0;
//	caseCounts=(char *)calloc(noCases,sizeof(char));
    
	
    controlCounts=new unsigned short[noControls];
    for(int i=0;i<noControls;i++)
        controlCounts[i]=0;
//	controlCounts=(char *)calloc(noControls,sizeof(char));
    
    pVal=0;
    significanceType='n';
}

Kmer::~Kmer()
{
//    delete [] caseCounts;
//    delete [] controlCounts;
}

void Kmer::freeMemory()
{
	delete [] caseCounts;
	delete [] controlCounts;
//	free(caseCounts);
//	free(controlCounts);
}


void Kmer::show()
{
    char kmerString[100];
    
    cout<<getKmer(kmer,kmerString,31)<<" ";
    
    for(int k=0;k<noCases;k++)
    {
        cout<<(int)caseCounts[k]<<" ";
        
    }
    for(int k=0;k<noControls;k++)
    {
        cout<<(int)controlCounts[k]<<" ";
        
    }
    cout<<pVal<<" "<<significanceType<<endl;
    
}


HashTable::HashTable()
{
    totalKmerCountsCase=new long long int[noCases];
    for(int i=0;i<noCases;i++)
        totalKmerCountsCase[i]=0;
    
    totalKmerCountsControl=new long long int[noControls];
    for(int i=0;i<noControls;i++)
        totalKmerCountsControl[i]=0;
	
	for(int i=0;i<HASH_TABLE_LENGTH;i++)
	{
		pthread_mutex_init(&hashTableBuckets_mutex[i],NULL);
		totalKmers_Bucket[i]=0;
	}
    
    totalKmers=0;
    totalTests=0;
}


void HashTable::insertKmer(long long int val, int count, int isCase, int sampleNo)
{
	unsigned long int index=getHash(val) % HASH_TABLE_LENGTH;
	int found=0;

	pthread_mutex_lock(&hashTableBuckets_mutex[index]);
    
	for(int i=0;i<kmers[index].size();i++)
	{
		if(val==kmers[index][i]->kmer)
        {
            if(isCase==1)
            {
                kmers[index][i]->caseCounts[sampleNo]=count;
            }
            else
            {
                kmers[index][i]->controlCounts[sampleNo]=count;
            }
			found=1;
			break;
		}
	}
	if(found==0)
	{
		Kmer* km=new Kmer(noCases,noControls);
		km->kmer=val;
		
        if(isCase==1)
        {
            km->caseCounts[sampleNo]=count;
        }
        else
        {
            km->controlCounts[sampleNo]=count;
        }
        
		kmers[index].push_back(km);
       
		totalKmers_Bucket[index]++;
        	

	}
	pthread_mutex_unlock(&hashTableBuckets_mutex[index]);
 
}

void HashTable::computePvalues()
{
    
    int kmerCountCase=0;
    int kmerCountControl=0;
    
    for(int i=0;i<HASH_TABLE_LENGTH;i++)
    {
        for(int j=0;j<kmers[i].size();j++)
        {
            kmerCountCase=0;
            kmerCountControl=0;
            
            
            for(int k=0;k<noCases;k++)
            {
                if((int)kmers[i][j]->caseCounts[k]>0)
                    kmerCountCase++;
                
            }
            for(int k=0;k<noControls;k++)
            {
                if((int)kmers[i][j]->controlCounts[k]>0)
                    kmerCountControl++;
                
            }
            if(kmerCountControl==noControls && kmerCountCase==noCases)
            {
                kmers[i][j]->pVal=1;
            }
            else
            {
                kmers[i][j]->pVal=getPvalue_chi(kmerCountCase,kmerCountControl , noCases-kmerCountCase, noControls-kmerCountControl);
                totalTests++;
                
                if(kmerCountCase/(double)noCases < (kmerCountCase+kmerCountControl)/(double)(noCases+noControls))
                {
                    kmers[i][j]->significanceType='a';
                    
                }
                else if(kmerCountCase/(double)noCases > (kmerCountCase+kmerCountControl)/(double)(noCases+noControls))
                {
                    kmers[i][j]->significanceType='p';
                    
                }
            }
        }
        
        
    }
    
}

void * likelihoodRatio_thread(void *threadid)
{
	long tid;
   	tid = (long)threadid;
	double meanCase, meanControl, mean;
    
	double likelihoodNull, likelihoodAlt, likelihoodRatio;


	long long int noKmersCases=0, noKmersControls=0;	

	int presentCount=0;
	double presentRatio=0;



	for(int k=0;k<noCases;k++)
       {
           noKmersCases+=ht->totalKmerCountsCase[k];
                
       }
	for(int k=0;k<noControls;k++)
       {
           noKmersControls+=ht->totalKmerCountsControl[k];
                
       }



    for(int i=tid;i<HASH_TABLE_LENGTH;i+=NUM_THREADS)
    {
        for(int j=0;j<ht->kmers[i].size();j++)
        {
            
            meanCase=0;
            meanControl=0;
            mean=0;
       
		presentCount=0;	 

     
            for(int k=0;k<noCases;k++)
            {
                meanCase+=ht->kmers[i][j]->caseCounts[k];
			if(ht->kmers[i][j]->caseCounts[k]>0)
			{
				presentCount++;
			}                
            }
            for(int k=0;k<noControls;k++)
            {
                meanControl+=ht->kmers[i][j]->controlCounts[k];
                if(ht->kmers[i][j]->controlCounts[k]>0)
			{
				presentCount++;
			}	
            }
            
            mean=(meanCase+meanControl)/(double)(noKmersCases+noKmersControls);
            presentRatio=presentCount/(double)(noCases+noControls);
           
		if(presentRatio>=0.01 && presentRatio<=0.99)
	     	{
		 	ht->kmers[i][j]->forPCA='y';
	     	}
		else
		{
			ht->kmers[i][j]->forPCA='n';
		}          
        	if(presentRatio>=0.05)
	     	{
		 	ht->kmers[i][j]->forPval='y';
	     	}
		else
		{
			ht->kmers[i][j]->forPval='n';
		}
            likelihoodNull=0;
            likelihoodAlt=0;

		likelihoodAlt+=(poissonProb(meanCase, meanCase));
		likelihoodAlt+=(poissonProb(meanControl, meanControl));

		likelihoodNull+=(poissonProb(meanCase, mean*noKmersCases));
		likelihoodNull+=(poissonProb(meanControl, mean*noKmersControls));
		

/*            
            for(int k=0;k<noCases;k++)
            {
                likelihoodAlt+=log(poissonProb(kmers[i][j]->caseCounts[k], meanCase));
                likelihoodNull+=log(poissonProb(kmers[i][j]->caseCounts[k], mean));
                
                
            }
            for(int k=0;k<noControls;k++)
            {
                likelihoodAlt+=log(poissonProb(kmers[i][j]->controlCounts[k], meanControl));
                likelihoodNull+=log(poissonProb(kmers[i][j]->controlCounts[k], mean));
                
            }
 */           
            likelihoodRatio=likelihoodAlt-likelihoodNull;
            
		if(likelihoodRatio<0)
		{
			likelihoodRatio=0;
		}		

            double pVal=alglib::chisquarecdistribution(1, 2*likelihoodRatio);
            
            ht->kmers[i][j]->pVal=pVal;
	
		meanCase=meanCase;
            
            	meanControl=meanControl*noKmersCases/noKmersControls;
 		
		ht->kmers[i][j]->meanCase=meanCase;
		ht->kmers[i][j]->meanControl=meanControl;

		if(meanCase>meanControl)
		{
			ht->kmers[i][j]->significanceType='p';
		}
		else if(meanCase<meanControl)
		{
			ht->kmers[i][j]->significanceType='a';

		}
		else
		{
			ht->kmers[i][j]->significanceType='n';

		}
        }
    }

	pthread_exit(NULL);


}

void HashTable::computeLikelihoodRatios()
{
    


	pthread_t threads[NUM_THREADS];
	int rc;
	long t;
	void *status;
	for(t=0; t<NUM_THREADS; t++)
	{
		  rc = pthread_create(&threads[t], NULL, likelihoodRatio_thread, (void *)t);
		  if (rc){
			 exit(-1);
		  }
	}
	
	for(t=0; t<NUM_THREADS; t++) 
	{
		rc = pthread_join(threads[t], &status);
		if (rc) 
	  	{
         		exit(-1);
      		}
      	}

}

FILE *caseFile;
FILE *controlFile;
double pValThreshold;


void * dump_thread(void *threadid)
{
	long tid;
   	tid = (long)threadid;

	char kmerString[100];
    
    for(int i=tid;i<HASH_TABLE_LENGTH;i+=NUM_THREADS)
    {
        for(int j=0;j<ht->kmers[i].size();j++)
        {

		// for eigenstrat
		if(ht->kmers[i][j]->forPCA=='y' && rand()/((double)(RAND_MAX)+1)<0.01)	//ht->kmers[i][j]->pVal<=pValThreshold && rand()/((double)(RAND_MAX)+1)<0.001)
		{
			pthread_mutex_lock(&eigenFile_mutex);


			fprintf(eigenSNPFile,"%s\t%d\t%lf\t%d\n",getKmer(ht->kmers[i][j]->kmer, kmerString, 31),1,0.0,0);

			for(int k=0;k<noCases;k++)
            		{
                		fprintf(eigenGenoFile,"%d\t",ht->kmers[i][j]->caseCounts[k]>0?1:0);
                
            		}
            		for(int k=0;k<noControls;k++)
            		{
               		fprintf(eigenGenoFile,"%d\t",ht->kmers[i][j]->controlCounts[k]>0?1:0);
                
            		}	
			fprintf(eigenGenoFile,"\n");
			pthread_mutex_unlock(&eigenFile_mutex);

		}



            if(ht->kmers[i][j]->pVal<=pValThreshold && ht->kmers[i][j]->forPval=='y')
            {
                if(ht->kmers[i][j]->significanceType=='p')
                {
              //      fprintf(caseFile,"%s\t%d\t%d\n",getKmer(kmers[i][j]->kmer, kmerString, 31),kmers[i][j]->caseCounts[0],kmers[i][j]->controlCounts[0]);

			pthread_mutex_lock(&caseOutFile_mutex);

                      fprintf(caseFile,"%s\t%d\t%d\t%e\t",getKmer(ht->kmers[i][j]->kmer, kmerString, 31),(int)ht->kmers[i][j]->meanCase,(int)ht->kmers[i][j]->meanControl,ht->kmers[i][j]->pVal);

			for(int k=0;k<noCases;k++)
            		{
                		fprintf(caseFile,"%d\t",ht->kmers[i][j]->caseCounts[k]);
                
            		}
            		for(int k=0;k<noControls;k++)
            		{
               		fprintf(caseFile,"%d\t",ht->kmers[i][j]->controlCounts[k]);
                
            		}	
			fprintf(caseFile,"\n");
			pthread_mutex_unlock(&caseOutFile_mutex);


                }
                else if(ht->kmers[i][j]->significanceType=='a')
                {
               //     fprintf(controlFile,"%s\t%d\t%d\n",getKmer(kmers[i][j]->kmer, kmerString, 31),kmers[i][j]->caseCounts[0],kmers[i][j]->controlCounts[0]);
	
			pthread_mutex_lock(&controlOutFile_mutex);
       	       fprintf(controlFile,"%s\t%d\t%d\t%e\t",getKmer(ht->kmers[i][j]->kmer, kmerString, 31),(int)ht->kmers[i][j]->meanCase,(int)ht->kmers[i][j]->meanControl,ht->kmers[i][j]->pVal);
			for(int k=0;k<noCases;k++)
            		{
                		fprintf(controlFile,"%d\t",ht->kmers[i][j]->caseCounts[k]);
                
            		}
            		for(int k=0;k<noControls;k++)
            		{
               		fprintf(controlFile,"%d\t",ht->kmers[i][j]->controlCounts[k]);
                
            		}	
			fprintf(controlFile,"\n");
			pthread_mutex_unlock(&controlOutFile_mutex);
       	       
                }

                
            }
        }
        
    }
    
    for(int i=tid;i<HASH_TABLE_LENGTH;i+=NUM_THREADS)
    {
        for(int j=0;j<ht->kmers[i].size();j++)
        {
		ht->kmers[i][j]->freeMemory();
		delete ht->kmers[i][j];
        }
        ht->kmers[i].clear();
    }


	pthread_exit(NULL);

}


void HashTable::dumpKmers(double sigLevel)
{
    caseFile=fopen("case_out_wo_bonf.kmerDiff","a");
    controlFile=fopen("control_out_wo_bonf.kmerDiff","a");

	eigenGenoFile=fopen("gwas_eigenstratX.geno","a");
    	eigenSNPFile=fopen("gwas_eigenstratX.snp","a");


	
	pValThreshold=sigLevel/(double)CUTOFF;    


	pthread_t threads[NUM_THREADS];
	int rc;
	long t;
	void *status;
	for(t=0; t<NUM_THREADS; t++)
	{
		  rc = pthread_create(&threads[t], NULL, dump_thread, (void *)t);
		  if (rc){
			 exit(-1);
		  }
	}
	
	for(t=0; t<NUM_THREADS; t++) 
	{
		rc = pthread_join(threads[t], &status);
		if (rc) 
	  	{
         		exit(-1);
      		}
      	}

	fclose(caseFile);
	fclose(controlFile);

    	fclose(eigenSNPFile);
	fclose(eigenGenoFile);


}

void HashTable::show()
{
/*    
    for(int i=0;i<HASH_TABLE_LENGTH;i++)
    {
        for(int j=0;j<kmers[i].size();j++)
        {
            kmers[i][j]->show();
        }
        
    }
*/
    cout<<totalKmers<<endl;
    cout<<totalTests<<endl;
}

Factorials::Factorials(int max)
{
    factorials=new double[max+1];
    double product=1;
    factorials[0]=1;
    for(int i=1;i<=max;i++)
    {
        product=product*i;
        factorials[i]=product;
    }
    
}

double Factorials::getFactorial(int number)
{
    return factorials[number];
    
}

		

void getKeyVal(char *s, KeyVal* kv)
{
	long long int	val=0;
	int countVal=0;
	int i=0;
	char ch;
	while(i<KMER_LENGTH)
	{
		
		ch=s[i++];
		
		val=val<<2|bases[ch];		
        
	}
	i++;
	kv->val=val;
	while(1)
	{
		ch=s[i];
		if(ch=='\0'||ch=='\n')
		{	
			break;
		}
		countVal=countVal*10+ch-'0';
		i++;
	        
	}
	kv->count=countVal;

}

FILE ** kmerFilesCases;
FILE ** kmerFilesControls;	
long long int *valsCases;
long long int *valsControls;
unsigned short int * countsCases;
unsigned short int * countsControls;

struct ThreadArg
{
	long long int valBar;
	int threadID;
};

void * readCases(void *threadid)
{
	ThreadArg *ta=(ThreadArg *)threadid;
	long long int valBar=ta->valBar;
   	int threadNo=ta->threadID;
	long long int val;
	int count;

	char *line= new char[MAX_REC_LEN];
    	int MAX_FILE_READ=MAX_REC_LEN/sizeof(line[0]);

	for(int i=threadNo;i<noCases;i+=NUM_THREADS/2)
    		{
			if(valsCases[i]!=-1 && valsCases[i]<valBar)
			{
				ht->insertKmer(valsCases[i], countsCases[i], 1, i);
				valsCases[i]=-1;
				countsCases[i]=-1;
			}
			if(valsCases[i]==-1)
			{
       		while(fscanf(kmerFilesCases[i],"%lld %d\n",&val,&count)!=EOF)
        		{
			
							
				if(val<valBar)
				{
 			           	ht->insertKmer(val, count, 1, i);
            			}
				else
				{
					valsCases[i]=val;
					countsCases[i]=count;
					break;
				}
        		}
 			}
    		}


	delete []line;
	pthread_exit(NULL);

}

void * readControls(void *threadid)
{
	ThreadArg *ta=(ThreadArg *)threadid;
	long long int valBar=ta->valBar;
   	int threadNo=ta->threadID;
	long long int val;
	int count;
	
	char *line= new char[MAX_REC_LEN];
    	int MAX_FILE_READ=MAX_REC_LEN/sizeof(line[0]);


	for(int i=threadNo;i<noControls;i+=NUM_THREADS/2)
    		{
			if(valsControls[i]!=-1 && valsControls[i]<valBar)
			{
				ht->insertKmer(valsControls[i], countsControls[i], 0, i);
				valsControls[i]=-1;
				countsControls[i]=-1;
			}
			if(valsControls[i]==-1)
			{
       		while(fscanf(kmerFilesControls[i],"%lld %d\n",&val,&count)!=EOF)
        		{

            			if(val<valBar)
				{
 			           	ht->insertKmer(val, count, 0, i);
            			}
				else
				{
					valsControls[i]=val;
					countsControls[i]=count;
					break;
				}
			
        		}
 			}
    		}

	delete []line;
	pthread_exit(NULL);

}

void printHelp()
{
    
	cout<<"hawk-"<<VERSION<<endl;
	cout<<"----------------"<<endl;
	cout<<endl;
	cout<<"hawk - hitting associations with k-mers"<<endl;
	cout<<"Usage:"<<endl;
	cout<<"hawk <noCases> <noControls>"<<endl;
	cout<<endl;
	cout<<"It requires counting k-mers with Jellyfish before running."<<endl;
	cout<<"Please see README for more details."<<endl;
	exit(1);
    
}


int main(int argc, const char * argv[])
{
	if(argc<2)
	{
		printHelp();
	}
	if(strcmp(argv[1],"--help")==0 || strcmp(argv[1],"-h")==0)
     		printHelp();	    

    noCases=atoi(argv[1]);
    noControls=atoi(argv[2]);
    
 
    ht=new HashTable();

	FILE *countsFile=fopen("case_total_kmers.txt","r");
	 
	for(int i=0;i<noCases;i++)
    	{
		fscanf(countsFile,"%lld\n",&ht->totalKmerCountsCase[i]);
	}
	countsFile=fopen("control_total_kmers.txt","r");
	
	for(int i=0;i<noControls;i++)
    	{
		fscanf(countsFile,"%lld\n",&ht->totalKmerCountsControl[i]);
	}


	FILE *caseFile=fopen("case_out_wo_bonf.kmerDiff","w");
    	FILE *controlFile=fopen("control_out_wo_bonf.kmerDiff","w");
	

	fclose(caseFile);
	fclose(controlFile);

	FILE *eigenGenoFile=fopen("gwas_eigenstratX.geno","w");
    	FILE *eigenSNPFile=fopen("gwas_eigenstratX.snp","w");

	fclose(eigenSNPFile);
	fclose(eigenGenoFile);



    
	kmerFilesCases=new FILE*[noCases];
	kmerFilesControls=new FILE*[noControls];	
    	char *kmerFilename;
    	kmerFilename=new char[5000];
    
	valsCases=new long long int[noCases];
	valsControls=new long long int[noControls];
	countsCases=new unsigned short int[noCases];
	countsControls=new unsigned short int[noControls];

	char *line= new char[MAX_REC_LEN];
	char *line2= new char[MAX_REC_LEN];
    	int MAX_FILE_READ=MAX_REC_LEN/sizeof(line[0]);



    
    char *temp;
    char kmerString[100];
    unsigned short count;
	
	long long int valMax=0x3FFFFFFFFFFFFFFF;
	long long int valBar=0;
	long long int val=0;
	long long int valInc=0x0010000000000000;

	FILE *sortedFile=fopen("case_sorted_files.txt","r");
	 
	for(int i=0;i<noCases;i++)
    	{
		fscanf(sortedFile,"%s\n",kmerFilename);
       //	sprintf(kmerFilename,"%d_case_kmers_sorted.txt",(i+1));
       	kmerFilesCases[i]=fopen(kmerFilename,"r");

		if(kmerFilesCases[i]==NULL)
		{
			cout<<kmerFilename<<" file doesn't exist"<<endl;
		}

		valsCases[i]=-1;
		countsCases[i]=-1;
		
	}
	sortedFile=fopen("control_sorted_files.txt","r");
	for(int i=0;i<noControls;i++)
    	{
		fscanf(sortedFile,"%s\n",kmerFilename);
 //      	sprintf(kmerFilename,"%d_control_kmers_sorted.txt",(i+1));
       	kmerFilesControls[i]=fopen(kmerFilename,"r");

		if(kmerFilesControls[i]==NULL)
		{
			cout<<kmerFilename<<" file doesn't exist"<<endl;
		}

		valsControls[i]=-1;
		countsControls[i]=-1;
		
	}

	
	ThreadArg *thArgsCase[NUM_THREADS/2];
	ThreadArg *thArgsControl[NUM_THREADS/2];

	while(valBar<valMax)
	{
		valBar+=valInc;


		pthread_t caseThreads[NUM_THREADS/2];
		pthread_t controlThreads[NUM_THREADS/2];
		int rc;
		void *status;
		for(int i=0;i<NUM_THREADS/2;i++)
		{
			thArgsCase[i]=new ThreadArg;
			thArgsCase[i]->valBar=valBar;
			thArgsCase[i]->threadID=i;
			rc = pthread_create(&caseThreads[i], NULL, readCases, (void *)thArgsCase[i]);
			if (rc)
			{
				 exit(-1);
			}
		}
		for(int i=0;i<NUM_THREADS/2;i++)
		{
			thArgsControl[i]=new ThreadArg;
			thArgsControl[i]->valBar=valBar;
			thArgsControl[i]->threadID=i;
			rc = pthread_create(&controlThreads[i], NULL, readControls, (void *)thArgsControl[i]);
			if (rc)
			{
				 exit(-1);
			}
		}
		
		for(int i=0;i<NUM_THREADS/2;i++)
		{
			rc = pthread_join(caseThreads[i], &status);
			delete thArgsCase[i];
			if (rc) 
	  		{
         			exit(-1);
      			}
		}
      		for(int i=0;i<NUM_THREADS/2;i++)
		{
			rc = pthread_join(controlThreads[i], &status);
			delete thArgsControl[i];
			if (rc) 
	  		{
         			exit(-1);
      			}
		}


		
	    	ht->computeLikelihoodRatios();
    
    	
    		ht->dumpKmers(SIG_LEVEL);
		
		cout<<valBar<<endl;
 
	}

	for(int i=0;i<HASH_TABLE_LENGTH;i++)
	{
		ht->totalKmers+=ht->totalKmers_Bucket[i];
		
	}

    	ht->show();




	caseFile=fopen("case_out_wo_bonf.kmerDiff","r");
    	controlFile=fopen("control_out_wo_bonf.kmerDiff","r");


    FILE *caseFileOut=fopen("case_out_w_bonf.kmerDiff","w");
    FILE *controlFileOut=fopen("control_out_w_bonf.kmerDiff","w");

	FILE *totalKmerFile=fopen("total_kmers.txt","w");
	fprintf(totalKmerFile,"%lld\n",ht->totalKmers);

	int count1, count2;
	double pVal;

    	 while(fgets(line, MAX_FILE_READ, caseFile)!=NULL)
        {
		strcpy(line2,line);


            	temp=strtok(line,"\t\n ");
            	strcpy(kmerString,temp);
            
            	temp=strtok(NULL,"\t\n ");
            	count1=atoi(temp);

		temp=strtok(NULL,"\t\n ");
            	count2=atoi(temp);
            			
		temp=strtok(NULL,"\t\n ");
            	pVal=atof(temp);

		if(pVal<SIG_LEVEL/ht->totalKmers)
		{
//			fprintf(caseFileOut,"%s\t%d\t%d\n",kmerString,count1,count2);
			fputs(line2,caseFileOut);

               
		}
	}


	while(fgets(line, MAX_FILE_READ, controlFile)!=NULL)
        {
		strcpy(line2,line);


            	temp=strtok(line,"\t\n ");
            	strcpy(kmerString,temp);
            
            	temp=strtok(NULL,"\t\n ");
            	count1=atoi(temp);

		temp=strtok(NULL,"\t\n ");
            	count2=atoi(temp);
            			
		temp=strtok(NULL,"\t\n ");
            	pVal=atof(temp);

		if(pVal<SIG_LEVEL/ht->totalKmers)
		{
		//	fprintf(controlFileOut,"%s\t%d\t%d\n",kmerString,count1,count2);
               	fputs(line2,controlFileOut);

		}
	}

    
    return 0;
}



