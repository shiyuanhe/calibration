#ifndef LINK_FUN_CLASS
#define LINK_FUN_CLASS

#include "RcppArmadillo.h"
using namespace arma;

class linkBase{
public:
    virtual ~linkBase(){
    }
    
    virtual vec theta_trans(const vec & theta) = 0;
    virtual vec theta_trans_deriv(const vec & theta) = 0;
    
    virtual vec theta_trans_deriv2O(const vec & theta) = 0;
    
    void set_bound(double lower_, double upper_){
        lowerB = lower_;
        upperB = upper_;
        rangeB = upperB - lowerB;
    }
    
protected:
    double lowerB, upperB, rangeB;
};


class linkIdentity: public linkBase{
public:
    vec theta_trans(const vec & theta){
        return theta;
    };
    vec theta_trans_deriv(const vec & theta){
        return ones(theta.n_elem);
    };
    
     vec theta_trans_deriv2O(const vec & theta){
        return zeros(theta.n_elem);
    }
    
    
};


class linkLogit: public linkBase{
public:
    vec theta_trans(const vec & theta){
        vec theta_new;
        theta_new = 1.0 + exp(-theta);
        theta_new = rangeB / theta_new + lowerB;
        return theta_new;
    };
    
    vec theta_trans_deriv(const vec & theta){
        vec theta_deriv_multi;
        theta_deriv_multi = 1.0 + exp(-theta);
        theta_deriv_multi = exp(-theta) / (theta_deriv_multi % theta_deriv_multi) * rangeB;
        return theta_deriv_multi;
    };
    
    vec theta_trans_deriv2O(const vec & theta){
        vec tmpExp, tmp1, result;
        tmpExp = exp(-theta);
        tmp1 = 1.0 + tmpExp;
        result = (tmpExp % tmpExp) / (tmp1 % tmp1 % tmp1) * (2.0 * rangeB);
        result -=  tmpExp / (tmp1 % tmp1)  * rangeB;
        return result;
    }
    
    
};



class expLink: public linkBase{
public:
    vec theta_trans(const vec & theta){
        return exp(theta);
    };
    
    vec theta_trans_deriv(const vec & theta){
        return exp(theta);
    };
    
    vec theta_trans_deriv2O(const vec & theta){
        return exp(theta);
    }
    
    
};

#endif

