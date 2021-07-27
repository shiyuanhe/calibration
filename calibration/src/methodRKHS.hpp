#ifndef CLASS_COREMETHOD
#define CLASS_COREMETHOD

#include "methodBase.hpp"

class calModel: public methodBase{
public:
    calModel(int c1, int c2, int c3, int numComp): methodBase(c1, c2, c3, 1){}

    void setLambda(double lambda_){
            thetaCoreVec[0]->set_lambda(lambda_);
    }
    
    
    // ojbective and gradient function
    double objFun(vec gamma);
    vec gradFun(vec gamma);
    mat predict_y(vec gamma, mat newX, bool getSigma);
    vec compute_sigma(const vec & gamma, const mat &newX);
    vec compute_sigma_version2(const vec & gamma, const mat &newX);
    mat fittedCovMat(vec gamma);
    double compute_GCV(vec gamma);
    vec compute_residualSigma(vec gamma);
private:
    // First order adjustment of sigma after mapping f(x)
    // X: sigma ---> F(x): |f'(x)|sigma
    vec adjust_sigma_1O(const vec & sigmaVec, const vec & gradVec);
    // Second order adjustment of sigma after mapping f(x)
    // X: sigma ---> F(x): |f'(x)|sigma
    vec adjust_sigma_2O(const vec & sigmaVec, const vec & gradVec,
                        const vec & hess);
    
};

#endif



