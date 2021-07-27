#ifndef THETACUBICV1_CLASS
#define THETACUBICV1_CLASS

#include "theta_RKHS.hpp"


class thetaRKSH_cubic_v1: public thetaRKSH{
public:
    thetaRKSH_cubic_v1(){
        set_dimNull(2);
    }
    

    double k2(const double &x){
        double result, k1;
        k1 = std::abs(x) - 0.5;
        result = 0.5 * (k1*k1 - 1.0/12.0);
        return result;
    }
    
    double k4(const double &x);
    
    mat computeKernel(const mat & dataX1, const mat & dataX2);
    
    mat computeNull(const mat & dataX1);
};


#endif
