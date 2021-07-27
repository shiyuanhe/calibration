#include "theta_fun.hpp"


vec thetaQuadratic::theta_value(const vec & gamma){
    vec result(sampleSize);
    result.fill(gamma(0));
    for(int i = 1; i <= dim_p; i++)
        result +=  gamma(i) * dataXScaled.col(i - 1) + 
            gamma(dim_p + i) * dataXScaled_Squared.col(i - 1);
    return result;
}

vec thetaQuadratic::theta_predict(const vec & gamma, const mat & newX){
    vec result(newX.n_rows);
    mat newXScaled, newXSquared;
    newXScaled = scaleData(newX);
    newXSquared = newXScaled % newXScaled;
    
    result.fill(gamma(0));
    for(int i = 1; i <= dim_p; i++)
        result +=  gamma(i) * newXScaled.col(i - 1) + 
            gamma(dim_p + i) * newXSquared.col(i - 1);
    
    return result;
}

mat thetaQuadratic::theta_deriv(const vec & gamma){
    mat result;
    result = join_vert(ones(sampleSize).t(), dataXScaled.t());
    result = join_vert(result, dataXScaled_Squared.t());
    return result;
}


mat thetaQuadratic::theta_deriv_newX(const vec & gamma, const mat & newX){
    mat result;
    mat newXScaled = scaleData(newX);
    mat newXSquared = newXScaled % newXScaled;
    result = join_vert(ones(newX.n_rows).t(), newXScaled.t());
    result = join_vert(result, newXSquared.t());
    return result;
}

void thetaQuadratic::set_data(const mat & dataX_, vec lowerB_, vec upperB_){
    thetaBase::set_data(dataX_, lowerB_, upperB_);
    dataXScaled_Squared = dataXScaled % dataXScaled;
    dim_p = dataXScaled.n_cols;
}
