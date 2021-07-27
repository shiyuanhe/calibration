#ifndef CLASS_BASEMETHOD
#define CLASS_BASEMETHOD

#include "RcppArmadillo.h"
#include "link_class.hpp"
#include "simuObj.hpp"
#include "simuObj_GP.hpp"

#include "theta_fun.hpp"
#include "theta_RKHS_cubic_v1.hpp"
#include "theta_RKHS_linear_v1.hpp"

#include <cmath>
#include <vector>
#include <memory>

using namespace std;
using namespace arma;
using namespace Rcpp;

class methodBase{
public:
    
    /**
    * c1: objective function
    * 0: y(theta, x) = theta
    *  1, 2, 3: simulation 1,2,3 
    *  0: GP
    */
    /**
    * c2: theta function
    * 1 Constant
    * 2 Exponential
    * 3 Qudratic
    * 4 RKHS Linear V1(Integral) V2(min max)
    * 5 RKHS Cubic V1(Integral) V2(min max)
    */
    /*
    * C3: Link Function 
    * 1: Identity link
    * 2: Logit Link
    */
    methodBase(int c1, int c2, int c3, int numComp_);
    virtual ~methodBase() {}
    
    void selectSimulator(int c1);
    void selectThetaType(int i, int c2);
    void selectLinkType(int i, int c3);
    int getNumComp() { return numComp;}
    
    // When C1 = 0, set data/beta for underlying GP
    void GP_setData(const vec & Y_, const mat X_){
        simulator->GP_setData(Y_, X_);
    }
    
    void GP_setBeta(const vec beta){
        simulator->GP_setBeta(beta);
    }
    
    void GP_setScaling(const vec scaling_){
        simulator->GP_setScaling(scaling_);
    }
    
    void GP_setNugget(const double sigmaSq){
        simulator->GP_setNugget(sigmaSq);
    }
    
    // Upper and lower bound for the link transformation
    virtual void set_linkBound(double lowerB_, double upperB_);
    
    // observed data, dataX has observed lower and upper bound
    virtual void setData(vec dataY_, mat dataX_, vec lowerB_, vec upperB_);
    
    // void setAuxData(mat aux_);
    // void setLambda(double lambda_){
    //     thetaCore->set_lambda(lambda_);
    // }
    // 
    // // ojbective and gradient function
    // double objFun(vec gamma);
    // vec gradFun(vec gamma);
    // mat predict_y(vec gamma, vec newX);
    
protected:
    int sampleSize; 
    vec dataY; // the response
    mat dataXReal; // data where the spline is applied on theta(x)
    // mat aux; // auxilary data for the obj function
    
    std::shared_ptr<simuobj> simulator;
    // std::shared_ptr<thetaBase> thetaCore;
    // std::shared_ptr<linkBase> thetaLink;
    size_t numComp; // the number of multivariate component
    std::vector<std::shared_ptr<thetaBase> > thetaCoreVec;
    std::vector<std::shared_ptr<linkBase> > thetaLinkVec;
    
};

#endif



