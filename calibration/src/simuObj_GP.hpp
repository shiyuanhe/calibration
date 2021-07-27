#ifndef SIMUOBJGP_CLASS
#define SIMUOBJGP_CLASS

#include "RcppArmadillo.h"
#include "simuObj.hpp"
#include "gaussianProcess.hpp"
using namespace arma;



// the GP simulation setup
class simuobj_GP:
    public simuobj{
public:
    
    void GP_setData(const vec & Y_, const mat X_){ internalGP.set_data(Y_, X_); }
    void GP_setBeta(const vec beta_){ internalGP.set_beta(beta_); }
    void GP_setScaling(vec scaling_){ internalGP.set_scaling(scaling_); }
    void GP_setNugget(double sigmaSq){ internalGP.set_nugget(sigmaSq); }
    vec yModel(const mat & theta, const mat & dataXReal){
        vec result;
        mat Xnew = join_horiz(dataXReal, theta);
        result = internalGP.predict_y(Xnew);
        return result;
    }
    
    mat yModelPartial(const mat & theta, const mat & dataXReal){
        int numSample, numComp, offSet, p;
        mat result;
        
        numSample = dataXReal.n_rows;
        offSet = dataXReal.n_cols;
        numComp = theta.n_cols;
        result = mat(numSample, numComp);
        mat Xnew = join_horiz(dataXReal, theta);
        
        for(p = 0; p < numComp; p++)
            result.col(p) = internalGP.emul_partial(Xnew, p + offSet);
        return result;
    }

    
    mat yModelPartial2O(const mat & theta, const mat & dataXReal){
        vec result;//, theta;
        mat Xnew = join_horiz(dataXReal, theta);
        result = internalGP.emul_partial2O(Xnew, Xnew.n_cols - 1);
        return result;
    }
    
    vec yModelPartial2OCross(const mat & theta, const mat & dataXReal,
                             int p, int q){
        vec result;//, theta;
        int offset = dataXReal.n_cols;
        mat Xnew = join_horiz(dataXReal, theta);
        result = internalGP.emul_partial2OCross(Xnew, p + offset, q + offset);
        return result;
    }
    
private:
    gaussianP internalGP;
};


#endif

