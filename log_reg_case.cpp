#include <cstring>
#include <pthread.h>
#include <semaphore.h>
#include <iostream>
#include <fstream>
#include <vector>

#include "lr.h"
#include "specialfunctions.h"

#define eps 1e-30

#define PCA_COUNT 10
//#define NULL_MODEL_FEATURE_COUNT 3
//#define ALT_MODEL_FEATURE_COUNT 4
#define MAX_LINE_LENGTH 512
#define CLASS_NAME_LENGTH 256
#define CLASS_NAME1 "Case"
#define CLASS_NAME2 "Control"
#define MALE "M"
#define FEMALE "F"
#define CHUNK_SIZE 10000
#define GENDER_INFO_FILE_NAME "gwas_info.txt"
// #define DEBUG


using namespace std;


// variables, semaphores and mutex used to sync and control thread life
pthread_mutex_t done_count_lock;
int done_count;
int thread_exit_signal;
sem_t all_start;
sem_t all_done;

ifstream con_file;
ifstream feature_z_file;
ifstream ind_file;
ifstream case_total_file;
ifstream control_total_file;
ifstream cov_file;
ifstream gender_info_file;

// These variables are global to ease the passing to multiple threads
int start_indx;
int num_of_thread;
int PC;
int cov_count;
string covfile;
string gender_info_file_name;
int read_row_count;

int NULL_MODEL_FEATURE_COUNT;
int ALT_MODEL_FEATURE_COUNT;

int mx_iter;
double learn_rate;

std::vector<std::vector<double> > Z;
std::vector<std::vector<double> > C;
std::vector<double> Y;
std::vector<unsigned long long int> totals;
std::vector<unsigned long long int> case_totals;
std::vector<unsigned long long int> control_totals;
std::vector<std::vector<unsigned long long int> > kmercounts;
std::vector<double> output;
std::vector<std::vector<double> > global_features_NULL;
std::vector<std::vector<double> > global_features_ALT;
std::vector<double> null_model;


struct thread_info
{
	int thread_no;
};


int open_file_connection();
int find_row_count();
void init_sync_primitve();
void * worker_thread_func(void *thread_no_as_ptr);

void printmodel(vector<double>model){
    for(size_t i = 0; i<model.size(); ++i){
        cout<<model[i]<<" ";
    }
    cout<<"\n";
}

int main(int argc,char **argv)
{
	num_of_thread = 1;
	PC = 2;
	covfile = "";
	cov_count = 0;
	gender_info_file_name = "";

	for(int i = 0; i<argc; ++i){
		if(strcmp(argv[i],"-t")==0){
			num_of_thread = atoi(argv[i+1]);
		}
		else if(strcmp(argv[i],"-c")==0){
			int ll = strlen(argv[i+1]);
			for(int j = 0; j<ll; ++j){
				covfile.push_back(argv[i+1][j]);
			}
		}
		else if(strcmp(argv[i],"-p")==0){
			PC = atoi(argv[i+1]);
		}
		else if(strcmp(argv[i],"-s")==0){
			gender_info_file_name = GENDER_INFO_FILE_NAME;
		}
	}

    if(open_file_connection()) {
        cout<<"Error in opening file"<<std::endl;
        return 0;
    }

    learn_rate = 0.1; 
    mx_iter = 25;

    unsigned int nrow = find_row_count();

	//Y.size() and nrow are equal
    Z = std::vector<std::vector<double> >(nrow,std::vector<double>(PCA_COUNT,0));
    Y = std::vector<double>(nrow);
    std::vector<double>RAWC;
    totals = std::vector<unsigned long long int>(nrow);
    #ifdef DEBUG
    std::vector<int> gender_info_ref(nrow);
    #endif
	for(unsigned int l=0;l<Y.size();l++)
    {
        char buf[MAX_LINE_LENGTH];
        char class_name[CLASS_NAME_LENGTH];
        for(int l1=0;l1<PCA_COUNT;l1++)
        {
            feature_z_file>>Z[l][l1];
        }

        /*
         * because each line's 3rd string (which is last) is of importance
         * just take three string input from the line and keep the last
         */
        ind_file.getline(buf,MAX_LINE_LENGTH-1);
        stringstream sstream(buf);
        // works only if class name is 3rd non alpha numeric seperated string of line\t'
		sstream>>class_name;
		#ifdef DEBUG
		std::cout<<class_name<<'\t';
		#endif
		sstream>>class_name;
		#ifdef DEBUG
		std::cout<<class_name<<'\t';
		if(strcmp(class_name,MALE)==0){
			gender_info_ref[l] = 1;
		}
		else if(strcmp(class_name,FEMALE)==0){
			gender_info_ref[l] = 0;
		}
		#endif
		sstream>>class_name;
		#ifdef DEBUG
		std::cout<<class_name<<std::endl;
		#endif
        if(strcmp(class_name,"Case")==0) {
            Y[l] = 1.0;
        }
        else {
            Y[l] = 0.0;
        }

        // as case_total_kmers.txt may have less entry than case_kmer_counts.txt+control_kmer_counts.txt 
		if(!case_total_file.eof()){
			unsigned long long int tmp;
			case_total_file>>tmp;
			if(!case_total_file.eof() && !case_total_file.fail()){
				case_totals.push_back(tmp);
			}
		}
		// as control_total_kmers.txt may have less entry than case_kmer_counts.txt+control_kmer_counts.txt
		if(!control_total_file.eof()){
			unsigned long long int tmp;
			control_total_file>>tmp;
			if(!control_total_file.eof() && !control_total_file.fail()){
				control_totals.push_back(tmp);
			}
		}
    }
    // construct totals
    for(unsigned int l=0;l<case_totals.size();l++)
	{
		totals[l]=case_totals[l];
	}
	for(unsigned int l=0;l<control_totals.size();l++)
	{
		totals[l+case_totals.size()]=control_totals[l];
	}
	//release the file connections
    feature_z_file.close();
    ind_file.close();
    case_total_file.close();
    control_total_file.close();

	// extracting sequence information
	// needed to sequence sex confounder info from gwas_info.txt and covariate info 
	// according to the sequence of samples in gwas_eigenstratX.ind
	// this would be probably cleaner if we used map<sample_id, class_id> from the begining
	std::vector<int> case_indx_info;
	std::vector<int> control_indx_info;
    char token[512];
	for(int line_no=0;line_no < nrow;line_no++){
		gender_info_file>>token;
		gender_info_file>>token;
		gender_info_file>>token;
		if(strcmp(token,CLASS_NAME1)==0){
			case_indx_info.push_back(line_no);
		}
		else if(strcmp(token,CLASS_NAME2)==0){
			control_indx_info.push_back(line_no);
		}
	}
	// rewinding gwas_info.txt file stream for future read
	gender_info_file.clear(); // needed before seekg if not c++11
    gender_info_file.seekg(0,ios::beg); 


    //reading covariate file...
    //like Z, i dunno how much PC is there
    if((int)covfile.size()>0){
    	double cc;
    	while(cov_file>>cc){
    		RAWC.push_back(cc);
    	}
    	int sz = (int)RAWC.size();
    	// first store them according to gwas_info.txt provided sample sequence
    	std::vector<std::vector<double> > C_tmp(nrow,std::vector<double>(sz/nrow,0));
    	int k = 0;
    	for(int i = 0; i<nrow; ++i){
    		for(int j = 0; j<(sz/nrow); ++j){
    			C_tmp[i][j] = RAWC[k];
    			k++;
    		}
    	}
    	cov_count = sz/nrow;
    	cov_file.close();
    	
    	// now store them in sequence according to gwas_eigenstratX.ind
		// first all case samples (stably sort) then the control samples (stably sort) information are filled
		// inside C 
		// we use earlier individually stored sequence of case sample index and sequence of control sample index
		C = std::vector<std::vector<double> >(nrow,std::vector<double>(sz/nrow,0));
		int C_fillup_indx = 0;
		for(int i = 0; i<case_indx_info.size(); ++i){
			int indx = case_indx_info[i];
			for(int j = 0; j<C_tmp[indx].size(); ++j){
				C[C_fillup_indx][j] = C_tmp[indx][j];
			}
			C_fillup_indx++;
		}
		for(int i = 0; i<control_indx_info.size(); ++i){
			int indx = control_indx_info[i];
			for(int j = 0; j<C_tmp[indx].size(); ++j){
				C[C_fillup_indx][j] = C_tmp[indx][j];
			}
			C_fillup_indx++;
		}
    }

	// similarly read gender covariate
    // reading gender covariate from gender info file of each read seq
    // assumption is same no. of line corresponds to same read seq. both in feature_z file
    // and gender info file
	// if only one seq. gender is unknown (not MALE or FEMALE) this info will not be used
    int unknown_gender = nrow;
    std::vector<int> gender_info;
    if((int)gender_info_file_name.size()>0){
		std::vector<int> gender_info_tmp(nrow);
		gender_info = std::vector<int>(nrow);
		
		char token[512];
		for(int line_no=0;line_no < nrow;line_no++){
			gender_info_file>>token;
			gender_info_file>>token;
			if(strcmp(token,MALE)==0){
				gender_info_tmp[line_no] = 1;
				unknown_gender--;
			}
			else if(strcmp(token,FEMALE)==0){
				gender_info_tmp[line_no] = 0;
				unknown_gender--;
			}

			gender_info_file>>token;
		}
		
		int gender_fillup_indx = 0;
		for(int i = 0; i<case_indx_info.size(); ++i){
			int indx = case_indx_info[i];
			gender_info[gender_fillup_indx] = gender_info_tmp[indx];
			gender_fillup_indx++;
		}
		for(int i = 0; i<control_indx_info.size(); ++i){
			int indx = control_indx_info[i];
			gender_info[gender_fillup_indx] = gender_info_tmp[indx];
			gender_fillup_indx++;
		}		
	}

    #ifdef DEBUG
    cout<<"gender infor file name "<<gender_info_file_name<<std::endl;
    cout<<"unknown gender count " << unknown_gender << std::endl;
    cout<<"nrow : "<<nrow<<std::endl;
    cout<<"Z"<<std::endl;
    for(int l=0;l<nrow;l++)
    {
        for(int l1=0;l1<PCA_COUNT;l1++)
        {
            cout<<Z[l][l1]<<' ';
        }
        cout<<std::endl;
    }

    cout<<"Y"<<std::endl;
    for(int l=0;l<nrow;l++)
    {
        cout<<Y[l]<<' ';
    }
    cout<<std::endl;

    cout<<"totals"<<std::endl;
    for(int l=0;l<nrow;l++)
    {
        cout<<totals[l]<<' ';
    }
    
    cout<<"covariate count " << cov_count <<std::endl;
    for(int i=0;i<C.size();i++)
    {
		for(int j = 0; j<C[i].size(); ++j){
    		std::cout<<C[i][j]<<' ';
    	}
    	std::cout<<std::endl;
    }
    cout<<std::endl;
    
    if(case_indx_info.size() + control_indx_info.size() == C.size()) {
		std::cout << "covariate sample count and case+control sample count match" << std::endl;
	}
    
    if(gender_info_ref.size() == gender_info.size()) {
		std::cout << "gender information sample count and case+control sample count match" << std::endl;
	}
    
    int correct = 1;
    for(int i=0;i<gender_info_ref.size();i++)
    {
		if(gender_info_ref[i]!=gender_info[i]){
			correct = 0;
		}
	}
	std::cout << "gender_info correctness " << correct << std::endl;
    
    #endif

	/*
	 * below matrix creation is done for fitting glm using glm function
	 * 4th column of matrix will be different for each sample (as per understanding)
	 */
	int chunk_size = CHUNK_SIZE;
	
	NULL_MODEL_FEATURE_COUNT = 1+PC+cov_count+1;
	// if gender_info read and every seq. gender is known it means gender info is read. 
	// So, one more feature
	if(gender_info_file_name.size()>0 && unknown_gender==0){
		NULL_MODEL_FEATURE_COUNT++;
	}
	ALT_MODEL_FEATURE_COUNT = 1+NULL_MODEL_FEATURE_COUNT;

	global_features_NULL = std::vector<std::vector<double> >(nrow,std::vector<double>(NULL_MODEL_FEATURE_COUNT));
	global_features_ALT = std::vector<std::vector<double> >(nrow,std::vector<double>(ALT_MODEL_FEATURE_COUNT));

	for(unsigned int l=0;l<nrow;l++)
	{
		global_features_NULL[l][0] = 1;
		global_features_ALT[l][0] = 1;
		for(unsigned int z = 0; z<PC; ++z){
			global_features_NULL[l][z+1] = Z[l][z];
			global_features_ALT[l][z+1] = Z[l][z];
		}
		for(unsigned int c = 0; c<cov_count;++c){
			global_features_NULL[l][1+PC+c] = C[l][c];
			global_features_ALT[l][1+PC+c] = C[l][c];
		}
		if(unknown_gender==0){
			global_features_NULL[l][1+PC+cov_count] = gender_info[l];
			global_features_ALT[l][1+PC+cov_count] = gender_info[l];
			global_features_NULL[l][1+PC+cov_count+1] = totals[l];
			global_features_ALT[l][1+PC+cov_count+1] = totals[l];
		}
		else{
			global_features_NULL[l][1+PC+cov_count] = totals[l];
			global_features_ALT[l][1+PC+cov_count] = totals[l];
		}
	}

	null_model = glm(global_features_NULL,Y,0.1,20);

	output = std::vector<double>(CHUNK_SIZE);
	int chunkread = 0;
    
    /*
     * First we will init sync primitive to create a signal for threads 
     * to start
     */
    init_sync_primitve();
    // Then create the thread(s)
    std::vector<pthread_t> thread_list(num_of_thread);
    for(int l=0;l<num_of_thread;l++)
    {
		thread_info *info = new thread_info;
		info->thread_no = l;
		pthread_create(&thread_list[l], NULL, worker_thread_func, (void *)info);
	}

    while(true)
    {
        char buf[MAX_LINE_LENGTH];
        kmercounts.clear();
        for(read_row_count=0;read_row_count<chunk_size;read_row_count++)
        {
			con_file>>buf>>buf>>buf>>buf;

			if(con_file.eof()) {
				break;
			}

			kmercounts.push_back(std::vector<unsigned long long int>(Y.size()));
			for(unsigned int l=0;l<Y.size();l++)
			{
				con_file>>kmercounts[read_row_count][l];
			}
		}

		//kmercounts er size protibar CHUNK_SIZE kore bartese
		//kmercounts er each row te 15 ta column

        #ifdef DEBUG_THREAD
        /*
         * loop to see extraction from con_file is done correctly
         */
        cout<<"Portion of kmercounts : "<<std::endl;
        for(int l=0;l<chunk_size;l++)
        {
			for(int l1=0;l1<Y.size();l1++)
			{
				cout<<kmercounts[l][l1]<<' ';
			}
			cout<<std::endl;
		}
        #endif
        //read_row_count is equal to CHUNK_SIZE
        //shudhu sesh bar ektu kom hoite pare
        /*
         * Below for loop can be done in parallel. The plan is to divide the  
         * available iteration to multiple thread. If one thread is complete 
         * it will increase a signal variable from sequence of signal and go to sleep. 
         * If all signals are marked than main thread will write the result in stdout then  
         * read more data and signal the worker thread to restart. Main thread will act 
         * as a watcher.
         * 
         * Every thread will do interleaved reading from kmercounts and write to 
         * particular loaction exclusive to thread. As it's write operation is not in 
         * same memory address for different thread no synchrnization needed
         */ 
        
        pthread_mutex_lock(&done_count_lock);
        /* 
         * done_count holds how many thread has completed their part
         * done_count is 0. So that, last thread to complete this iteration 
         * can read the done_count and signal main thread appropriately
		 */
		done_count = 0;
		pthread_mutex_unlock(&done_count_lock);

		sem_init(&all_done, 0, 0);
		// sem_init(&all_start, 0, num_of_thread-1);
		for(int l=0;l<num_of_thread;l++) sem_post(&all_start);
        
        sem_wait(&all_done);
		// all_done semphore is posted, means all outputs are prepared. 
		// So, dump them to stdout
		
		for(int l=0;l<read_row_count;l++)
		{
			cout<<output[l]<<endl;
		}
		//cout<<"write"<<endl;
		cout.flush();
		//cout<<"flush"<<endl;
		
        if(read_row_count < chunk_size) {
			thread_exit_signal = 1;
			for(int l=0;l<num_of_thread;l++) sem_post(&all_start);
			break;
		}
		start_indx += read_row_count;

    }

    return 0;
}


int open_file_connection()
{
    con_file.open("case_out_w_bonf_top.kmerDiff");
    feature_z_file.open("pcs.evec");
    ind_file.open("gwas_eigenstratX.ind");
    case_total_file.open("case_total_kmers.txt");
    control_total_file.open("control_total_kmers.txt");
    if((int)covfile.size()>0){
    	char cvv[200];
    	for(int i = 0; i<(int)covfile.size();++i){
    		cvv[i] = covfile[i];
    		cvv[i+1] = '\0';
    	}
    	cov_file.open(cvv);
    	if(!cov_file){
    		cout<<covfile<<" not found";
    		return 1;
    	}
    }
    
    if((int)gender_info_file_name.size()>0){
		char gifn[200];
		for(int i = 0; i<(int)gender_info_file_name.size();i++){
			gifn[i] = gender_info_file_name[i];
		}
		gifn[gender_info_file_name.size()] = '\0';
		gender_info_file.open(gifn);
		if(!gender_info_file){
			cout<<gender_info_file_name<<" not found";
			return 1;
		}
	}

    if(!con_file) {
        cout<<"case_out_w_bonf_top.kmerdiff not found";
        return 1;
    }
    if(!feature_z_file) {
        cout<<"pcs.evec not found";
        return 1;
    }
    if(!ind_file) {
        cout<<"gwas_eigenstratX.ind not found";
        return 1;
    }
    if(!case_total_file) {
        cout<<"case_total_kmers.txt not found";
        return 1;
    }
    if(!control_total_file) {
        cout<<"control_total_kmers.txt not found";
        return 1;
    }
    
    return 0;
}


int find_row_count()
{
    int nrow = 0;
    char buf[MAX_LINE_LENGTH];

    while(!case_total_file.eof())
    {
        case_total_file.getline(buf,MAX_LINE_LENGTH-1);
        if(!case_total_file.eof()){
			nrow++;
		}
    }
    while(!control_total_file.eof())
    {
        control_total_file.getline(buf,MAX_LINE_LENGTH-1);
        if(!control_total_file.eof()){
			nrow++;
		}
    }
    case_total_file.clear(); // needed before seekg if not c++11
    case_total_file.seekg(0,ios::beg);
	control_total_file.clear(); // needed before seekg if not c++11
    control_total_file.seekg(0,ios::beg);

	return nrow;
}


void init_sync_primitve()
{
	sem_init(&all_start, 0, 0);
	sem_init(&all_done, 0, 0);
	pthread_mutex_init(&done_count_lock, 0);
	thread_exit_signal = 0;
}


void * worker_thread_func(void *arg)
{
	int thread_no = ((thread_info *)arg)->thread_no;
	int interleave = num_of_thread;
	
	std::vector<double> counts(Y.size());
	std::vector<std::vector<double> > thread_local_features_ALT(global_features_ALT);
	
	while(true)
	{
		sem_wait(&all_start);
		if(thread_exit_signal==1) {
			break;
		}
		
		for(int l=thread_no;l<read_row_count;l+=interleave)
		{
			for(unsigned int l1=0;l1<Y.size();l1++)
			{
				counts[l1] = kmercounts[l][l1]/(double)totals[l1];
			}
			
			//create the fourth column of matrix
			for(unsigned int l1=0;l1<Y.size();l1++)
			{
				thread_local_features_ALT[l1][ALT_MODEL_FEATURE_COUNT-1] = counts[l1];
			}
			
			std::vector<double> alt_model = glm(thread_local_features_ALT,Y,learn_rate,mx_iter);
			
			double alt_likelihood = 1.0;
			for(int dat = 0; dat < thread_local_features_ALT.size(); ++dat){
				std::vector<double> data(thread_local_features_ALT[0].size());
				for(int j = 0; j<thread_local_features_ALT[0].size(); ++j){
					data[j] = thread_local_features_ALT[dat][j];
				}
				double p = predict(alt_model, data);
				
				if( Y[dat] == 1){
					alt_likelihood = alt_likelihood*p;
				}
				else{
					alt_likelihood = alt_likelihood*(1.0-p); 
				}
			}

			double null_likelihood = 1.0;
			for(int dat = 0; dat < global_features_NULL.size(); ++dat){
				std::vector<double> data(global_features_NULL[0].size());
				for(int j = 0; j<global_features_NULL[0].size(); ++j){
					data[j] = global_features_NULL[dat][j];
				}
				double p = predict(null_model, data);

				if( ((int)Y[dat]) == 1){
					null_likelihood = null_likelihood*p;
				}
				else{
					null_likelihood = null_likelihood*(1.0-p); 
				}
			}

			if((null_likelihood)==0 && (alt_likelihood)==0.0) {
				null_likelihood = 0.001;
				alt_likelihood = 1.0;
			}

			double likelihood_ratio = null_likelihood/alt_likelihood;
			double log_likelihood_ratio = -2.0*(log(likelihood_ratio));

			if(fabs(log_likelihood_ratio)<eps || log_likelihood_ratio<0.0){
				log_likelihood_ratio = 0.0;
			}

			output[l] = alglib::chisquarecdistribution(1, log_likelihood_ratio);
		}
		
		pthread_mutex_lock(&done_count_lock);
		done_count++;
		if(done_count==num_of_thread) {
			sem_post(&all_done);
		}
		pthread_mutex_unlock(&done_count_lock);
	}
	
	// delete ((thread_info *)arg);
}
