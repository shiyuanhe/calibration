#ifndef GPROCESS_CLASS
#define GPROCESS_CLASS

#include "RcppArmadillo.h"
#include "squaredExp.hpp"
#include <cmath>
using namespace std;
using namespace arma;


// For all related dataX input, each obs is stored in the rows.
// In this class, the data is transposed to store
// in the columns.
class gaussianP{
public:
    // Set scaling, the default values are all one
    void set_scaling(vec scaling_){
        if(hasData) kernelObj.set_scaling(scaling_);
    }
    
    double cv_obj(vec beta);
    vec cv_grad(vec beta);
    
    vec predict_y(mat Xnew);
    mat predict_yCov(mat Xnew);
    
    vec emul_partial(mat Xnew, int c1);

    vec emul_partial2O(mat Xnew, int c1);
    vec emul_partial2OCross(mat Xnew, int c1, int c2); 
    void set_beta(vec beta);
    
    void set_data(vec y_, mat X_){
        y = y_;
        X = X_.t();
        nSample = X.n_cols;
        hasData = true;
        kernelObj.set_scaling(ones(X.n_rows));
    }
    
    void set_nugget(double sigmaSq){
        nugget = sigmaSq;
    }
    
    gaussianP();
private:

    int nSample;
    squaredExp kernelObj;
    // The matrix X stores obs in columns.
    bool hasData;
    mat X, K, Kinv;
    vec y, alpha, kernelBeta;
    double nugget;
    
    double cv_grad_j(vec beta, int j);
    
};

#endif
