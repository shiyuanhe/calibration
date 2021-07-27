#include "theta_RKHS_cubic_v1.hpp"


double thetaRKSH_cubic_v1::k4(const double &x){
    double result, k1, k1Sq;
    k1 = std::abs(x) - 0.5;
    k1Sq = k1 * k1;
    result = k1Sq * k1Sq - 0.5 * k1Sq + 7.0/240.0;
    result /= 24.0;
    return result;
}


mat thetaRKSH_cubic_v1::computeKernel(const mat & dataX1, const mat & dataX2){
    int ns1, ns2;
    ns1 = dataX1.n_rows;
    ns2 = dataX2.n_rows;
    
    mat kernelMat = mat(ns1, ns2, fill::zeros);
    double diff;
    for(int i = 0; i< ns1; i++)
        for(int  j = 0; j < ns2; j++){
            diff = dataX1(i) - dataX2(j);
            kernelMat(i,j) = k2(dataX1(i)) * k2(dataX2(j)) - k4(diff);
        }
        return kernelMat;
}

mat thetaRKSH_cubic_v1::computeNull(const mat & dataX1){
    int dataN = dataX1.n_rows;
    return join_horiz(ones(dataN), dataX1);
}
