#ifndef THETA_CLASS
#define THETA_CLASS

#include "RcppArmadillo.h"
using namespace arma;

class thetaBase{
public:
    virtual ~thetaBase(){
        
    }
    
    thetaBase(){
        isRKHS = false;
    }
    
    bool get_isRKHS(){ 
        return isRKHS; 
    }
    
    virtual vec theta_value(const vec & gamma) = 0;
    
    // each column for one observation
    // each row for one variable. 
    virtual mat theta_deriv(const vec & gamma) = 0;
    virtual mat theta_deriv_2O(const vec & gamma, const vec coef){
        return mat(gamma.n_elem, gamma.n_elem, fill::zeros);
    }
    
    virtual vec theta_predict(const vec & gamma, const mat & newX) = 0;
    virtual mat theta_deriv_newX(const vec & gamma, const mat & newX) = 0;
    
    
    // default: no penalty
    virtual double penalty(const vec & gamma){
        return 0.0;
    }
    
    virtual vec penalty_deriv(const vec & gamma){
        return zeros(gamma.n_elem);
    }
    
    virtual mat penalty_deriv2O(const vec & gamma){
        return mat(gamma.n_elem, gamma.n_elem, fill::zeros);
    }
    
    
    void set_lambda(double lambda_){
        lambda = lambda_;
    }
    
    virtual void set_data(const mat & dataX_, vec lowerB_, vec upperB_){
        upperB = upperB_;
        lowerB = lowerB_;
        dataXScaled = scaleData(dataX_);
        sampleSize = dataXScaled.n_rows;
    }
    
    virtual double compute_residual_sigma(const vec & eta, const vec & Y,
                                          const vec & ym,  const vec & ym1D){
        return 0.1;
    }
    virtual double compute_GCV(const vec & eta, const vec & Y,
                               const vec & ym,  const vec & ym1D){
        return 0.0;
    }
    virtual vec RKHS_theta_sigma(const mat & newX, const vec & ym1D, double sigmaErr){
        return (-1.0 * ones(newX.n_rows));
    }
    
    mat scaleData(const mat & dataX_){
        int i, nr, nc;
        nr = dataX_.n_rows;
        nc = dataX_.n_cols;
        mat tmp(nr, nc);
        for(i = 0; i < nc; i++){
            tmp.col(i) = dataX_.col(i);
            tmp.col(i) -= lowerB(i);
            tmp.col(i) /= (upperB(i) - lowerB(i));
        }
        return tmp;
    }
    
protected:
    mat dataXScaled;
    int sampleSize;
    double lambda;
    vec upperB, lowerB;
    bool isRKHS;
    
    
};


class thetaConst: public thetaBase{
public:
    vec theta_value(const vec & gamma){
        vec result;
        result = ones(sampleSize) * gamma(0);
        return result;
    }
    
    vec theta_predict(const vec & gamma, const mat & newX){
        vec result;
        result = ones(newX.n_rows) * gamma(0);
        return result;
        
    }
    
    mat theta_deriv(const vec & gamma){
        return ones(sampleSize).t();
    }
    
    mat theta_deriv_newX(const vec & gamma, const mat & newX){
        return ones(newX.n_rows).t();
    }
    
};


// gamma(0) * exp(-gamma(1)  x_1 - gamma(2) x_2 - ... - gamma(p) x_p)
class thetaExponential: public thetaBase{
public:
    vec theta_value(const vec & gamma);
    
    vec theta_predict(const vec & gamma, const mat & newX);
    
    mat theta_deriv(const vec & gamma);
    
    mat theta_deriv_newX(const vec & gamma, const mat & newX);
    
    mat theta_deriv_2O(const vec & gamma, const vec coef);
    
    void set_data(const mat & dataX_, vec lowerB_, vec upperB_);
    
private:
    int dim_p;
};


// gamma(0) + gamma(1) x_1 + ... + gamma(p) x_p
//   + gamma(p+1) x_1^2 + ... + gamma(2p) x_p^2
class thetaQuadratic: public thetaBase{
public:
    vec theta_value(const vec & gamma);
    
    vec theta_predict(const vec & gamma, const mat & newX);
    
    mat theta_deriv(const vec & gamma);
    
    mat theta_deriv_newX(const vec & gamma, const mat & newX);
    
    void set_data(const mat & dataX_, vec lowerB_, vec upperB_);
    
private:
    mat dataXScaled_Squared;
    int dim_p;
};

#endif
