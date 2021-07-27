#ifndef THETALINEARV1_CLASS
#define THETALINEARV1_CLASS

#include "theta_RKHS.hpp"

class thetaRKSH_linear_v1: public thetaRKSH{
public:
    thetaRKSH_linear_v1(){
        set_dimNull(1);
    }
    
private:

    inline double k2(const double &x){
        double result, k1;
        k1 = std::abs(x) - 0.5;
        result = 0.5 * (k1 * k1 - 1.0/12.0);
        return result;
    }
    
    inline double k1(const double &x){
        return std::abs(x) - 0.5;
    }
    
    mat computeKernel(const mat & dataX1, const mat & dataX2){
        int ns1, ns2;
        ns1 = dataX1.n_rows;
        ns2 = dataX2.n_rows;
        
        mat kernelMat = mat(ns1, ns2, fill::zeros);
        double  diff;
        for(int i = 0; i< ns1; i++)
            for(int  j = 0; j < ns2; j++){
                diff = dataX1(i) - dataX2(j);
                kernelMat(i,j) = k1(dataX1(i)) * k1(dataX2(j)) +
                    k2(diff);
            }
            
            return kernelMat;
    }
    
    mat computeNull(const mat & dataX1){
        int dataN = dataX1.n_rows;
        return ones(dataN);
    }
    
    
};

#endif
