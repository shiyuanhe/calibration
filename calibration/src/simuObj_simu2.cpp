#include "simuObj.hpp"
// Objective model for simulation 2


// model
vec simuobj_simu2::yModel(const mat & theta, const mat & dataXReal){
    vec result; 
    
    vec tmp;
    tmp =  0.5 *(dataXReal - 2.0) % (dataXReal - 2.0) + 0.5;
    tmp = theta / tmp;
    result = cos(2.0 * dataXReal) % sin(0.5 * dataXReal);
    double a = 0*PI, b = 0*PI, c = 3.0;
    result %= sin(a*tmp+0.5*PI) % cos(b*tmp) % exp(c*(tmp - 1.0));
    return result;
}

// model: gradient with respect to theta
mat simuobj_simu2::yModelPartial(const mat & theta, const mat & dataXReal){
    int nSample = dataXReal.n_rows;
    vec result, baseV, addV;
    vec tmp1, tmp2;
    double a = 0*PI, b = 0*PI, c = 3.0;
    
    tmp1 =  0.5 *(dataXReal - 2.0) % (dataXReal - 2.0) + 0.5;
    tmp1 = 1.0 / tmp1;
    tmp2 = theta % tmp1;
    baseV = cos(2.0 * dataXReal) % sin(0.5 * dataXReal);
    
    result = zeros(nSample);
    addV = a*tmp1 % cos(a*tmp2+0.5*PI) % cos(b*tmp2) % exp(c*tmp2-c);
    result += baseV % addV;
    
    addV =  -b*tmp1 % sin(a*tmp2+0.5*PI) % sin(b*tmp2) % exp(c*tmp2-c);
    result += baseV % addV;
    
    addV = c*tmp1 % sin(a*tmp2+0.5*PI) % cos(b*tmp2) % exp(c*tmp2-c);
    result += baseV % addV;
    
    return result;
}

// model: second order (2O) gradient with respect to theta
mat simuobj_simu2::yModelPartial2O(const mat & theta, const mat & dataXReal){
    int nSample = dataXReal.n_rows;
    vec result, baseV, addV;
    vec tmp1, tmp2, tmp1Sq;
    double a = 0*PI, b = 0*PI, c = 3.0;
    
    tmp1 =  0.5 * (dataXReal - 2.0) % (dataXReal - 2.0) + 0.5;
    tmp1 = 1.0 / tmp1;
    tmp1Sq = tmp1 % tmp1;
    tmp2 = theta % tmp1;
    baseV = cos(2.0 * dataXReal) % sin(0.5 * dataXReal);
    
    result = zeros(nSample);
    // addV = (PI/2.0) * tmp1 % cos((PI/2.0) * tmp2) % cos((2.0*PI)*tmp2) % exp(10.0 * tmp2 - 10.0);
    // result += baseV % addV;
    addV = (-a*a) * tmp1Sq % sin(a*tmp2+0.5*PI) % cos(b*tmp2) % exp(c*tmp2-c);
    result += baseV % addV;
    addV = (-a*b) * tmp1Sq % cos(a*tmp2+0.5*PI) % sin(b*tmp2) % exp(c*tmp2-c);
    result += baseV % addV;
    addV = (a*c) * tmp1Sq % cos(a*tmp2+0.5*PI) % cos(b*tmp2) % exp(c*tmp2-c);
    result += baseV % addV;
    
    // addV = (-2.0 * PI) * tmp1 % sin((PI/2.0) * tmp2) % sin((2.0*PI)*tmp2) % exp(10.0 * tmp2 - 10.0);
    // result += baseV % addV;
    addV =  (-b * a) * tmp1Sq % cos(a*tmp2+0.5*PI) % sin(b*tmp2) % exp(c*tmp2-c);
    result += baseV % addV;
    addV =  (-b* b) * tmp1Sq % sin(a*tmp2+0.5*PI) % cos(b*tmp2) % exp(c*tmp2-c);
    result += baseV % addV;
    addV =  (-b * c) * tmp1Sq % sin(a*tmp2+0.5*PI) % sin(b*tmp2) % exp(c*tmp2-c);
    result += baseV % addV;
    
    // addV = 10.0 * tmp1 % sin((PI/2.0) * tmp2) % cos((2.0*PI)*tmp2) % exp(10.0* tmp2 - 10.0);
    // result += baseV % addV;
    addV = (c*a) * tmp1Sq % cos(a*tmp2+0.5*PI) % cos(b*tmp2) % exp(c*tmp2-c);
    result += baseV % addV;
    addV = (-c*b) * tmp1Sq % sin(a*tmp2+0.5*PI) % sin(b*tmp2) % exp(c*tmp2-c);
    result += baseV % addV;
    addV = (c*c) * tmp1Sq % sin(a*tmp2+0.5*PI) % cos(b*tmp2) % exp(c*tmp2-c);
    result += baseV % addV;
    
    return result;
}
