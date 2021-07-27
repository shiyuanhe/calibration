#ifndef SQUARED_EXP_HEADER
#define SQUARED_EXP_HEADER

#include "RcppArmadillo.h"
#include <cmath>
using namespace std;
using namespace arma;


// All data sample in matrix columns

class squaredExp{
public:
    squaredExp(){
        beta = vec(2, fill::zeros);
        theta = exp(beta);
    }
    
    void set_beta(vec beta_){
        beta = beta_;
        theta = exp(beta);
    }
    
    void set_scaling(vec scaling_){
        scalingParam = scaling_;
    }
    
    //k = 0 covariance exp(beta[0]) * exp( - distSq / exp(beta[1]) )
    //k = 1 partial deriv wrt beta1, the value is the same as k = 0
    //k = 2 partial deriv wrt beta2, 
    mat kernelVar(const mat & X, const int k = 0);
    
    mat kernelCov(const mat & X1, const mat & X2);
    
    // The first order partial direvative w.r.t X1(cj,i)
    mat kernelCov_partial_cj(const mat & X1, const mat & X2, int cj);
    
    // The second order partial direvative w.r.t X1(cj,i)
    mat kernelCov_partial2O_cj(const mat & X1, const mat & X2, int cj);
    
    // The second order partial direvative w.r.t X1(cj,i)
    mat kernelCov_partial2OCross(const mat & X1, const mat & X2, 
                                 int c1, int c2);
    
    double compute_distSq(const vec & X1colI, const vec & X2colJ);
    
    double compute_kij(const vec & X1colI, const vec & X2colJ);
    
private:
    vec scalingParam;
    vec theta, beta;
};


#endif
