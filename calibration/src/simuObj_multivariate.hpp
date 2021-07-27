#ifndef SIMUOBJMULTI_CLASS
#define SIMUOBJMULTI_CLASS

#include "simuObj.hpp"


// the third simulation setup
class simuobj_multi1:
    public simuobj{
public:
    vec yModel(const mat & theta, const mat & dataXReal){
        vec result; 
        result = theta.col(0) % dataXReal;
        result += theta.col(1) % dataXReal % dataXReal;
        return result;
    }
    
    // The first order partial derivatives
    // Each row for one observation
    // Each column for one theta
    mat yModelPartial(const mat & theta, const mat & dataXReal){
        mat result(dataXReal.n_rows, 2);
        result.col(0) = dataXReal;
        result.col(1) = dataXReal % dataXReal;
        return result;
    }
    

    vec yModelPartial2OCross(const mat & theta, 
                             const mat & dataXReal,  
                             int i, int j){
        return arma::zeros<vec>(dataXReal.n_rows);
    }
    
};



// the fourth simulation setup
class simuobj_multi2:
    public simuobj{
public:
    vec yModel(const mat & theta, const mat & dataXReal){
        vec result(dataXReal.n_rows); 
        
        size_t i;
        for(i = 0; i < dataXReal.n_rows; i++){
            result(i) = theta(i,0) * pow(dataXReal(i,0), theta(i,1));
        }
        return result;
    }
    mat yModelPartial(const mat & theta, const mat & dataXReal){
        mat result(dataXReal.n_rows, 2);
        size_t i;
        for(i = 0; i < dataXReal.n_rows; i++){
            result(i,0) = pow(dataXReal(i,0), theta(i,1));
            result(i,1) = result(i,0) * theta(i,0) * log(dataXReal(i,0));
        }
        return result;
    }
    
    
    vec yModelPartial2OCross(const mat & theta, 
                             const mat & dataXReal, 
                             int p, int q){
        vec result(dataXReal.n_rows, fill::zeros); 
        size_t i;
        if(p==0 && q ==0){
            result = arma::zeros(dataXReal.n_rows, 1);
        }
        
        if(p==1 && q ==1){
            for(i = 0; i < dataXReal.n_rows; i++){
                result(i) = theta(i,0) * pow(dataXReal(i,0), theta(i,1))
                * std::log(dataXReal(i,0)) * std::log(dataXReal(i,0));
            }
        }
        
        if(p!=q){
            for(i = 0; i < dataXReal.n_rows; i++){
                result(i) = pow(dataXReal(i,0), theta(i,1)) *
                    std::log(dataXReal(i,0));
            }
        }
        
        
        return result;
    }
    
};




#endif
