#include <bits/stdc++.h>
#include "util.hpp"
#include "matrix_methods.h"

#define D_matrix std::vector<vector<double> >

using namespace std;

bool debug = true;

std::ostream & operator<<(std::ostream &stream,std::vector<double> v) 
{
	for(int l = 0;l<v.size();l++)
	{
		stream<<v[l]<<' ';
	}
	return stream;
}

double sigmoid(double x) {
    double e = 2.718281828;
    return 1.0 / (1.0 + pow(e, -x));
}

D_matrix _multiply(D_matrix &m1, D_matrix &m2){
	D_matrix ans;
	for(int i = 0; i<m1.size(); ++i){
		std::vector<double>line(m2[0].size(),0);
		ans.push_back(line);
	}
	if(m1[0].size()!=m2.size()){
		cout<<"cannot _multiply\n";
		return ans;
	}
	for(size_t i = 0; i<m1.size(); ++i){
		for(size_t j = 0; j<m2[0].size(); ++j){
			for(size_t k = 0; k<m1[0].size(); ++k){
				ans[i][j] = ans[i][j] + m1[i][k]*m2[k][j];
			}
		}
	}
	return ans;
}

double linear_predictor(std::vector<double>&model, std::vector<double>&data){
	double s = 0.0;
	for(int i = 0; i<model.size(); ++i)s = s + model[i]*data[i];

	return s;
}

double predict(std::vector<double>&model, std::vector<double>&data){
	double s = 0.0;
	for(int i = 0; i<model.size(); ++i)s = s + model[i]*data[i];

	return sigmoid(s);
}

D_matrix transpose(D_matrix &mm){
	D_matrix m(mm[0].size(),std::vector<double>(mm.size(),-1));
	for(int i = 0;i<mm.size(); ++i){
		for(int j = 0;j<mm[0].size();++j){
			m[j][i] = mm[i][j];
		}
	}
	return m;
}

void printD(D_matrix mm){
	for(size_t i = 0; i<mm.size(); ++i){
		for(size_t j = 0; j<mm[0].size(); ++j){
			cout<<mm[i][j]<<" ";
		}
		cout<<"\n";
	}
	return;
}

std::vector<double> glm_newton_raphson(std::vector<std::vector<double> >& x, std::vector<double>& y, double _gamma, int mx_iters, bool &is_singularity_error, bool &is_nan_error, double &ret_error, int &exit_iteration_no) {
    srand(time(0)); double epsilon = 1e-6;
    double gamma = _gamma;
    double error = 0.0;
    int max_iters = mx_iters;
    int iter = 0;

    std::vector<double> weight_old(x[0].size());
    for (size_t i=0; i<weight_old.size(); ++i) {
        double mxx = -10000000000.0;
        for(size_t j = 0; j<x.size(); ++j){
        	mxx = max(mxx,x[j][i]);
        }
        weight_old[i] = 1.0/mxx;
    }
    double prev_error = 1e18;
    D_matrix A = x; D_matrix A_T = transpose(A);
    D_matrix A_minus_y(A.size(),std::vector<double>(1,0));
    std::vector<double>B_proxy(A.size(),0);
    D_matrix A_T_B(A_T.size(),std::vector<double>(B_proxy.size(),0));

	while(true){
		double error = 0.0;
		for(size_t i = 0; i<A.size(); ++i){
			double z_i = 0.0;
			for(size_t j = 0; j<A[0].size(); ++j){
				z_i = z_i + A[i][j]*weight_old[j];
			}
			double alph_i = sigmoid(z_i);
			error = error + (y[i]-alph_i)*(y[i]-alph_i);
			B_proxy[i] = alph_i*(1.0-alph_i); 
			A_minus_y[i][0] = alph_i-y[i]; 
		}

		error /= x.size();
		ret_error = error;
		if (fabs(error-prev_error) < epsilon) {
			break;
        }

        prev_error = error;

		for(size_t i = 0; i<A_T_B.size(); ++i){
			for(size_t j = 0; j<A_T_B[0].size(); ++j){
				A_T_B[i][j] = A_T[i][j]*B_proxy[j];
			}
		}

		D_matrix toinv = _multiply(A_T_B,A);

		bool sing = false;
		bool nan = false;

		D_matrix hinv = inverse(toinv,(int)toinv.size(),sing,nan);
		if(sing || nan){
			is_singularity_error = sing;
			is_nan_error = nan;
			ret_error = prev_error;
			return weight_old;
		}

		D_matrix gradient = _multiply(A_T,A_minus_y);
		D_matrix grad_mul_hinv = _multiply(hinv,gradient);

		for(size_t j = 0; j<weight_old.size(); ++j){
			weight_old[j] = weight_old[j] - gamma*grad_mul_hinv[j][0];
		}

		iter += 1;
		exit_iteration_no = iter;
        if (iter >= max_iters) {
			break;
        }

        prev_error = error;
		ret_error = prev_error;
	}

	return weight_old;
}


/*
 * Logistic regression model optimization using IRLS algorithm
 * code is written following
 * http://114.70.193.188/BML/slides_freda/lec7.pdf
 * and
 * https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/glm.R
 * used function interface can be reduced to smaller interface
 */
std::vector<double> glm_irls(std::vector<std::vector<double> >& x, std::vector<double>& y, double _gamma, int mx_iters, bool &is_singularity_error, bool &is_nan_error, double &ret_error, int &exit_iteration_no) {
	double epsilon = 1e-6;
    int max_iters = mx_iters;
    int iter = 0;

	// initializing needed data structure early
	// maybe helpful in reducing allocation time cost

	// weight vector to return
    std::vector<double> weight(x[0].size(),1);
    // weight matrix of shape (feature count,1) considering offset coefficient as feauture also
    // this one will be deduced and assigned at each iteration
    D_matrix w(x[0].size(),std::vector<double>(1,1));
    // to store the good observation/sample from given x at each iteration
    D_matrix X(x.size(),std::vector<double>(x[0].size(),0));
    // to store the linear predictor \sum{x*w} at each iteration
    D_matrix eta(x.size(),std::vector<double>(1,0));
    // to store sigmoid(eta) at each iteration
    D_matrix mu(x.size(),std::vector<double>(1,0));
    // to store mu*(1-mu) of the good observation at each iteration
    // it was supposed to be diagonal matrix of (# of good observation, # of good observation)
    // made as (#of input observation, 1) to store only the diagonal element
    D_matrix S(x.size(),std::vector<double>(1,0));
    // to store fisher score(!)
    // z = eta + (y-mu)/(mu*(1-mu) + 1e-305) of the good observations at each iteration
    D_matrix z(x.size(),std::vector<double>(1,0));
    // to store S*X, resultant shape (# of good observation, # of feature)
    // considering offset coefficient as a feature also
    D_matrix S_X(x.size(),std::vector<double>(x[0].size(),0));
    // to store S*X, resultant shape (# of good observation, 1) at each iteration
    D_matrix S_z(x.size(),std::vector<double>(1,0));

    // in the iteration loop, weight will be deduced at end
    // loop will begin with assumption of mu and eta
    // this initialization is done following R function of  binomial$initialize
    // can be seen using following R code
    // > f = binomial()
	// > f = binomial(logit)
	// > f$initialize
	// basically mu is 1/4 if y is 0 else 3/2
	for(size_t i=0;i<mu.size();i++)
	{
		mu[i][0] = (y[i] + 0.5)/2;
		eta[i][0] = log(mu[i][0]/(1-mu[i][0]));
	}

	// initial value of square error which will be returned
    double prev_error = 1e18;
    // this loop may be broken by four condition
    // iteration > max_iteration
    // change in error < epsilon(1e-6)
    // matrix inversion is not possible
    // no good observation is not found
	while(true){
		// clearing vectors for new iteration
		X.clear();
		S.clear();
		S_X.clear();
		z.clear();
		S_z.clear();

		double error = 0.0;
		// good observation means observation with mu*(1-mu) or variance(!) > 0
		// 1e-305 is taken as epsilon
		int num_of_good_observations = 0;
		for(size_t i = 0; i<x.size(); ++i){
			double g_i = mu[i][0]*(1.0-mu[i][0]);
			// taking only the good observations
			if (g_i > 1e-305) {
				num_of_good_observations++;

				X.push_back(x[i]);

				z.push_back(std::vector<double>(1, eta[i][0] + (y[i]-mu[i][0])/(g_i+1e-305)));
				S.push_back(std::vector<double>(1,g_i));
			}

			error = error + (y[i]-mu[i][0])*(y[i]-mu[i][0]);
		}
		// if no good observation is found model fitting is exited
		// following https://github.com/SurajGupta/r-source/blob/a28e609e72ed7c47f6ddfbb86c85279a0750f0b7/src/library/stats/R/glm.R#L231
		if (num_of_good_observations == 0) {
			std::cerr << "No good observation is found" << std::endl;
			// weight got at previous iteration which is deduced from some good observation will be returned
			// this may be wrong because weight which produce no good observation is suspicious
			break;
		}

		// checking average error
		// only normal way of model fitting finishing
		error /= x.size();
		ret_error = error;
		if (fabs(error-prev_error) < epsilon) {
			break;
        }

        prev_error = error;

		// transpose for calculating hessian = X_T * S * X
		D_matrix X_T = transpose(X);
        // using simple loop instead of _multiply because S is supposed to be a diagonal matrix
        // reducing number of multiplication and needed memory
		for(size_t i=0;i<S.size();i++)
		{
			S_X.push_back(std::vector<double>(X[0].size(),0));
			for(size_t j=0;j<X[0].size();j++)
			{
				S_X[i][j] = S[i][0]*X[i][j];
			}
		}
		D_matrix hessian = _multiply(X_T, S_X);

		// to extract information that if hessian inversing cannot be done
		bool sing = false;
		bool nan = false;

		D_matrix hessian_inv = inverse(hessian,(int)hessian.size(),sing,nan);

		if(sing || nan){
			is_singularity_error = sing;
			is_nan_error = nan;

			ret_error = prev_error;
			break;
		}

		// to calculate X_T * S * z
		// using simple loop instead of _multiply because S is supposed to be a diagonal matrix
		// reducing number of multiplication
		for(size_t i=0;i<z.size();i++)
		{
			S_z.push_back(std::vector<double>(1,S[i][0]*z[i][0]));
		}
		D_matrix X_T_S_z = _multiply(X_T,S_z);
		// according to http://114.70.193.188/BML/slides_freda/lec7.pdf
		// updated w = (X_T * S * X)^-1 * X_T * S * z
		w = _multiply(hessian_inv,X_T_S_z);

		iter += 1;
		exit_iteration_no = iter;
        if (iter >= max_iters) {
           break;
        }

        prev_error = error;
		ret_error = prev_error;

		// assigning newly got w to variable which will be returned
		for(size_t i=0;i<weight.size();i++)
		{
			weight[i]=w[i][0];
		}

		// calculating eta and mu for next iteration using updated w
		for(size_t i = 0; i<x.size(); ++i){
			eta[i][0] = 0;
			for(size_t j = 0; j<x[0].size(); ++j){
				eta[i][0] += x[i][j]*w[j][0];
			}
			mu[i][0] = sigmoid(eta[i][0]);
		}
	}

	return weight;
}

std::vector<double> linear_regression(std::vector<std::vector<double> >& x, std::vector<double>& y) {

	std::vector< double > weight(x[0].size());
	for(int i = 0; i<(int)weight.size(); ++i){
		weight[i] = 0.0;
	}

    D_matrix X = x; D_matrix X_T = transpose(X);
    D_matrix _y = from_vector_to_D(y);

    D_matrix A = _multiply(X_T,X);
   	bool sing = false;
   	bool nan = false;
   	D_matrix A_inv = inverse(A,(int)A.size(),sing,nan);

   	if(sing){
		return weight;
   	}

   	D_matrix B = _multiply(A_inv,X_T);
   	D_matrix W = _multiply(B,_y);


    for(int i = 0; i<(int)W.size(); ++i){
    	weight[i] = W[i][0];
    }

	return weight;
}



