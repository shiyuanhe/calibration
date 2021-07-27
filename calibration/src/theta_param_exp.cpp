#include "theta_fun.hpp"

vec thetaExponential::theta_value(const vec & gamma){
    vec result;
    result = exp(- dataXScaled* gamma(span(1, dim_p))) * gamma(0);
    return result;
}

vec thetaExponential::theta_predict(const vec & gamma, const mat & newX){
    vec result;
    mat newXScaled;
    newXScaled = scaleData(newX);
    result = exp(- newXScaled* gamma(span(1, dim_p))) * gamma(0);
    return result;
}

mat thetaExponential::theta_deriv(const vec & gamma){
    mat row1, row2(sampleSize, dim_p);
    mat result;
    row1 = -gamma(0) * exp(-dataXScaled* gamma(span(1, dim_p)));
    for(int i = 0; i<dim_p; i++)
        row2.col(i) = row1 % dataXScaled.col(i);
    row1 /= -gamma(0);
    result = join_vert(row1.t(), row2.t());
    return result;
}

mat thetaExponential::theta_deriv_newX(const vec & gamma, const mat & newX){
    mat row1, row2(newX.n_rows, dim_p);
    mat result;
    mat newXScaled = scaleData(newX);
    
    row1 = -gamma(0) * exp(-newXScaled* gamma(span(1, dim_p)));
    for(int i = 0; i<dim_p; i++)
        row2.col(i) = row1 % newXScaled.col(i);
    row1 /= -gamma(0);
    result = join_vert(row1.t(), row2.t());
    return result;
}


mat thetaExponential::theta_deriv_2O(const vec & gamma, const vec coef){
    mat result(dim_p + 1, dim_p + 1, fill::zeros), tmp;
    tmp = result;
    double eValue;
    int c,i,j;
    for(c = 0; c < coef.n_elem; c++){
        eValue = - as_scalar(dataXScaled.row(c) * gamma(span(1, dim_p)));
        eValue = exp(eValue);
        tmp(0,0) = 0.0;
        for(i = 0; i < dim_p; i++){
            tmp(0,i + 1) = -eValue * dataXScaled(c,i);
            tmp(i + 1,0) = tmp(0, i + 1);
        }
        for(i = 0; i < dim_p; i++){
            for(j = i; j < dim_p; j++){
                tmp(i+1, j+1) = gamma(0) * eValue * 
                    dataXScaled(c,i) * dataXScaled(c,j);
                tmp(j+1, i+1) = tmp(i+1, j+1);
            }
        }
        result += tmp * coef(c);
    }
    
    return result;
}

void thetaExponential::set_data(const mat & dataX_, vec lowerB_, vec upperB_){
    thetaBase::set_data(dataX_, lowerB_, upperB_);
    dim_p = dataXScaled.n_cols;
}
