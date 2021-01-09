#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <cstring>
#include <vector>

using namespace std;
#define D_matrix std::vector<vector<double> >

D_matrix from_vector_to_D(std::vector<double>x){
    D_matrix dummy;
    for(int i = 0; i<(int)x.size(); ++i){
        std::vector<double>line;
        dummy.push_back(line);
        dummy[i].push_back(x[i]);
    }
    return dummy;
}


void printMatrix(D_matrix mm){
    for(size_t i = 0; i < mm.size(); ++i){
        for(size_t j = 0; j < mm[0].size(); ++j){
            cerr << setw(10) << mm[i][j] << "\t";
        }
        cerr << endl;
    }
    cerr << endl;
    return;
}

D_matrix initNewMatrix(int r, int c, double val){

    D_matrix newMat;
    for(int i = 0; i < r; i++){
        std::vector<double> emptyRow(c, val);
        newMat.push_back(emptyRow);
    }

    return newMat;
}

D_matrix initUnitMatrix(int r, int c){

    D_matrix newMat;
    for(int i = 0; i < r; i++){
        std::vector<double> emptyRow;
        for(int j = 0; j < c; j++){
            if(i == j){
                emptyRow.push_back(1);
            } else {
                emptyRow.push_back(0);
            }
        }
        newMat.push_back(emptyRow);
    }

    return newMat;
}

D_matrix initMatWithRandom(int r, int c) {
    
    D_matrix newMat;
    int factor1 = 7, factor2 = 3;
    for (int i = 0; i < r; i++) {
        vector<double> row;
        for (int j = 0; j < c; j++) {
            double data = (((i + 1) * factor1) + (j + 1) * factor2) / 4.2;
            row.push_back(data);
            factor1 += (rand() % 13);
            factor2 += (rand() % 7);
        }
        newMat.push_back(row);
    }

    return newMat;
}

D_matrix multiply(D_matrix m1, D_matrix m2){
    D_matrix ans;
    for(int i = 0; i<m1.size(); ++i){
        std::vector<double>line(m2[0].size(),0);
        ans.push_back(line);
    }
    if(m1[0].size()!=m2.size()){
        cout<<"cannot multiply\n";
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

// Doolittle algorithm
void luDecomposition(D_matrix matOriginal, int n, D_matrix &matLower, D_matrix &matUpper) { 
    
    matLower = initNewMatrix(n, n, 0);
    matUpper = initNewMatrix(n, n, 0);
    
    // Decomposing matrix into Upper and Lower 
    // triangular matrix 
    for (int i = 0; i < n; i++) { 
        // Upper Triangular 
        for (int k = i; k < n; k++) { 
            // Summation of L(i, j) * U(j, k) 
            double sum = 0; 
            for (int j = 0; j < i; j++) {
                sum += (matLower[i][j] * matUpper[j][k]); 
            }
            // Evaluating U(i, k) 
            matUpper[i][k] = matOriginal[i][k] - sum; 
        } 

        // Lower Triangular 
        for (int k = i; k < n; k++) { 
            if (i == k) 
                matLower[i][i] = 1; // Diagonal as 1 
            else { 

                // Summation of L(k, j) * U(j, i) 
                double sum = 0; 
                for (int j = 0; j < i; j++){
                    sum += (matLower[k][j] * matUpper[j][i]); 
                }

                // Evaluating L(k, i) 
                matLower[k][i] = (matOriginal[k][i] - sum) / matUpper[i][i]; 
            } 
        } 
    } 
} 

D_matrix inverse(D_matrix matOriginal, int n, bool &singular, bool &nan){

    D_matrix matLower, matUpper;
    luDecomposition(matOriginal, n, matLower, matUpper);

    D_matrix matInverse = initNewMatrix(n, n, 0);
    double det = 1;
     
    for(int inverse_col = 0; inverse_col < n; inverse_col++){

        vector<double> b;
        for(int j = 0; j < n; j++){
            if(inverse_col == j){
                b.push_back(1);
            } else {
                b.push_back(0);
            }
        }
        vector<double> y(n, 0);

        det *= matLower[0][0];
        // Forward substitution. Solve: Ly=b
        y[0] = b[0];
        for(int row = 1; row < n; row++){
            double sum = 0;
            for(int col = 0; col < n; col++){
                sum += matLower[row][col] * y[col];
            }
            y[row] = b[row] - sum;

            det *= matLower[row][row];
        }

        vector<double> x(n, 0);
        // Backward substitution. Solve: Ux=y
        x[n-1] = y[n-1] / matUpper[n - 1][n - 1];
        det *= matUpper[n - 1][n - 1];
        for(int row = n - 2; row > -1; row--){
            double sum = 0;
            for(int col = row + 1; col < n; col++){
                sum += matUpper[row][col] * x[col];
            }
            x[row] = (y[row] - sum) / matUpper[row][row];

            det *= matUpper[row][row];
        }

        for(int j = 0; j < n; j++){
            matInverse[j][inverse_col] = x[j];
        }
    }
	if(det == 0){
        singular = true;
    } else if(std::isnan(det)){
		nan = true;
	}
    else {
        singular = false;
        nan = false;
    }

    return matInverse;
}

