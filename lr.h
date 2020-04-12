#include <bits/stdc++.h>
/*#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>*/
#include "util.hpp"
#include "matrix_methods.h"
//#define BNU boost::numeric::ublas
#define D_matrix std::vector<vector<double> >
//namespace ublas = boost::numeric::ublas;

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

// template<class T>
// boost::numeric::ublas::matrix<T> gjinverse(const boost::numeric::ublas::matrix<T> &m, bool &singular){
//      using namespace boost::numeric::ublas;
//      const int size = m.size1();
//      // Cannot invert if non-square matrix or 0x0 matrix.
//      // Report it as singular in these cases, and return 
//      // a 0x0 matrix.
//      if (size != m.size2() || size == 0)
//      {
//          singular = true;
//          matrix<T> A(0,0);
//          return A;
//      }
//      // Handle 1x1 matrix edge case as general purpose 
//      // inverter below requires 2x2 to function properly.
//      if (size == 1)
//      {
//          matrix<T> A(1, 1);
//          if (m(0,0) == 0.0)
//          {
//              singular = true;
//              return A;
//          }
//          singular = false;
//          A(0,0) = 1/m(0,0);
//          return A;
//      }
//      // Create an augmented matrix A to invert. Assign the
//      // matrix to be inverted to the left hand side and an
//      // identity matrix to the right hand side.
//      matrix<T> A(size, 2*size);
//      matrix_range<matrix<T> > Aleft(A, range(0, size),range(0, size));
//      Aleft = m;
//      matrix_range<matrix<T> > Aright(A, range(0, size), range(size, 2*size));
//      Aright = identity_matrix<T>(size);
//      // Swap rows to eliminate zero diagonal elements.
//      for (int k = 0; k < size; k++)
//      {
//          if ( A(k,k) == 0 ) // XXX: test for "small" instead
//          {
//              // Find a row(l) to swap with row(k)
//              int l = -1;
//              for (int i = k+1; i < size; i++) 
//              {
//                  if ( A(i,k) != 0 )
//                  {
//                      l = i; 
//                      break;
//                  }
//              }
//              // Swap the rows if found
//              if ( l < 0 ) 
//              {
//                  //std::cerr << "Error:" <<  __FUNCTION__ << ":"
//                    //        << "Input matrix is singular, because cannot find"
//                      //      << " a row to swap while eliminating zero-diagonal.";
//                  singular = true;
//                  return Aleft;
//              }
//              else 
//              {
//                  matrix_row<matrix<T> > rowk(A, k);
//                  matrix_row<matrix<T> > rowl(A, l);
//                  rowk.swap(rowl);
//             #if defined(DEBUG) || !defined(NDEBUG)
//                  //std::cerr << __FUNCTION__ << ":"
//                    //        << "Swapped row " << k << " with row " << l 
//                      //      << ":" << A << "\n";
//             #endif
//              }
//          }
//      }
//      // Doing partial pivot
//      for (int k = 0; k < size; k++)
//      {
//          // normalize the current row
//          for (int j = k+1; j < 2*size; j++)
//              A(k,j) /= A(k,k);
//          A(k,k) = 1;
//          // normalize other rows
//          for (int i = 0; i < size; i++)
//          {
//              if ( i != k )  // other rows  // FIX: PROBLEM HERE
//              {
//                  if ( A(i,k) != 0 )
//                  {
//                      for (int j = k+1; j < 2*size; j++)
//                          A(i,j) -= A(k,j) * A(i,k);
//                      A(i,k) = 0;
//                  }
//              }
//          }
//         #if defined(DEBUG) || !defined(NDEBUG)
//          //std::cerr << __FUNCTION__ << ":"
//            //        << "GJ row " << k << " : " << A << "\n";
//         #endif
//      }
//      singular = false;
//      return Aright;
// }

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

std::vector<double> glm(std::vector<std::vector<double> >& x, std::vector<double>& y, double _gamma=0.1, int mx_iters=1000) {
    srand(time(0)); double epsilon = 0.000001;
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
		
        /*BNU::matrix<double>hessian(toinv.size(),toinv.size());
		for(size_t i = 0; i<toinv.size(); ++i){
			for(size_t j = 0; j<toinv.size(); ++j){
				hessian(i,j) = toinv[i][j];
			}
		}*/

		bool sing = false;
		//BNU::matrix<double> hessinv = gjinverse(hessian,sing);
		
		
		D_matrix hinv = inverse(toinv,(int)toinv.size(),sing);
		
		if(sing){
			return weight_old;
		}

		/*for(size_t i = 0; i<hinv.size(); ++i){
			for(size_t j = 0; j<hinv.size(); ++j){
				hinv[i][j] = hessinv(i,j);
			}
		}*/

		D_matrix gradient = _multiply(A_T,A_minus_y);
		D_matrix grad_mul_hinv = _multiply(hinv,gradient);

		for(size_t j = 0; j<weight_old.size(); ++j){
			weight_old[j] = weight_old[j] - gamma*grad_mul_hinv[j][0];
		}

		iter += 1;
        if (iter >= max_iters) {
            break;
        }
        
        prev_error = error;

	}

	return weight_old;
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
   	D_matrix A_inv = inverse(A,(int)A.size(),sing);

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



