#ifndef CLASS_COREMETHOD_MULTIVARIATE
#define CLASS_COREMETHOD_MULTIVARIATE

#include "methodBase.hpp"
#include "simuObj_multivariate.hpp"

class calModelMulti: public methodBase{
public:
    calModelMulti(int c1, int c2, int c3, int numComp);

    // ojbective and gradient function
    double objFun(mat gamma);
    mat gradFun(mat gamma);
    mat predict_y(mat gamma, mat newX, bool getSigma);
    mat sigmaSq_y(mat gamma, mat newX, double lambda);
    mat sigmaSq_y_simu(mat, mat , mat, double );
    
    double compute_GCV(mat gamma, double lambda);
    double compute_sigmaSq(mat gamma, double lambda);
    
    void setLambda(double lambda_){
        size_t i;
        for(i = 0; i < numComp; i++)
            thetaCoreVec[i]->set_lambda(lambda_);
    }
    
private:
    mat predict_yDeriv(mat gamma, mat newX);
    
    // compute the numerator and denominator
    // for the GCV and Sigma Squared
    vec compute_IminusA(mat gamma, double lambda);
    
    mat computeTheta(const mat & gamma);
    vec computeYHat(const mat & gamma);
    mat computeYHat_1O(const mat & gamma);
    
    // information matrix for the parametric model

    mat RKHS_sigmaY(mat gamma, mat newX, double lambda);
    vec Para_sigmaY(mat gamma, mat newX, double lambda);
    mat parametric_info(const mat & gamma);
};

#endif



