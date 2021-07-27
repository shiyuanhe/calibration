#ifndef SIMUOBJ_CLASS
#define SIMUOBJ_CLASS

#include "RcppArmadillo.h"
using namespace arma;


// computer model, objective and first order derivative
class simuobj{
protected:
    //vec dataXReal;
public:
    // virtual void set_data(const vec & dataXReal_){
    //     dataXReal = dataXReal_;
    // }
    virtual ~simuobj(){}
    virtual void GP_setData(const vec & Y_, const mat X_){}
    virtual void GP_setBeta(const vec beta){}
    virtual void GP_setScaling(vec scaling_){}
    virtual void GP_setNugget(double sigmaSq){}
    
    virtual vec yModel(const mat & theta,  const mat & dataXReal) = 0;
    // partial y(theta(x), x) against theta
    virtual mat yModelPartial(const mat & theta, const mat & dataXReal) = 0;
    virtual mat yModelPartial2O(const mat & theta, const mat & dataXReal){
        return arma::zeros(1,1);
    }
    virtual vec yModelPartial2OCross(const mat & theta, 
                                     const mat & dataXReal,
                                     int p, int q){
        return zeros<vec>(1);
    }
};


// y(theta(x), x) = theta(x)
class simuobj_null:
    public simuobj{
public:
    vec yModel(const mat & theta, const mat & dataXReal){
        return theta;
    }
    mat yModelPartial(const mat & theta, const mat & dataXReal){
        return arma::ones(theta.n_elem);
    }
    mat yModelPartial2O(const mat & theta, const mat & dataXReal){
        return arma::zeros(theta.n_elem);
    }
    
};


// the first simulation setup
class simuobj_simu1:
    public simuobj{
public:
    vec yModel(const mat & theta, const mat & dataXReal){
        vec result, coreF; 
        coreF = exp(dataXReal / 10.0) % cos(dataXReal);
        result = 0.5 * exp(dataXReal/5.0);
        result =  coreF % result / theta;
        return result;
    }
    mat yModelPartial(const mat & theta, const mat & dataXReal){
        vec result, coreF;
        coreF = exp(dataXReal / 10.0) % cos(dataXReal);
        result = 0.5 * exp(dataXReal/5.0);
        result =  - coreF % result / (theta % theta);
        return result;
    }
    
    mat yModelPartial2O(const mat & theta, const mat & dataXReal){
        vec coreF, result;
        coreF = exp(dataXReal / 10.0) % cos(dataXReal);
        result =  0.5 * exp(dataXReal/5.0);
        result =  2.0 * coreF % result / (theta % theta % theta);
        return result;
    }
    
};


// the second simulation setup
class simuobj_simu2:
    public simuobj{
public:
    vec yModel(const mat & theta, const mat & dataXReal);
    mat yModelPartial(const mat & theta, const mat & dataXReal);
    mat yModelPartial2O(const mat & theta, const mat & dataXReal);
    
};




#endif
