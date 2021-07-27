#ifndef THETA_RKHS_BASE_CLASS
#define THETA_RKHS_BASE_CLASS

#include "theta_fun.hpp"
#include <memory>
using namespace arma;

class thetaRKSH: public thetaBase{
public:
    virtual ~thetaRKSH(){}
    thetaRKSH(){
        isRKHS = true;
    }
    void set_data(const mat & dataX_,  vec lowerB_, 
                  vec upperB_);
    
    vec theta_value(const vec & gamma);
    
    vec theta_predict(const vec & gamma, const mat & newX);
    
    
    mat theta_deriv(const vec & gamma){
        return KmatLarge;
    }
    
    mat theta_deriv_newX(const vec & gamma, const mat & newX){
        mat result;
        mat newXScaled = scaleData(newX);
        mat KCov = computeKernel(newXScaled, dataXScaled);
        
        if(dimNull > 0){
            mat SCov = computeNull(newXScaled);
            result = join_horiz(SCov, KCov);
            result = result.t();
        }
        
        return result;
    }
    
    
    
    double penalty(const vec & gamma){
        double res;
        vec beta = gamma(span(dimNull, dimTotal - 1));
        res = lambda * as_scalar(beta.t() * Kmat * beta);
        return res;
    }
    
    vec penalty_deriv(const vec & gamma){
        vec res(dimTotal, fill::zeros);
        res(span(dimNull, dimTotal - 1)) += 
            Kmat *gamma(span(dimNull, dimTotal - 1));
        res = (2*lambda) * res;
        return res;
    }
    
    
    mat penalty_deriv2O(const vec & gamma){
        mat res(dimTotal, dimTotal, fill::zeros);
        res(span(dimNull, dimTotal - 1), span(dimNull, dimTotal - 1)) = Kmat;
        res = (2*lambda) * res;
        return res;
    }
    
    
    //vec RKHS_sigma(const vec & newX,const vec & ym1D, double sigmaErr);
    double compute_residual_sigma(const vec & eta, const vec & Y,
                                  const vec & ym,  const vec & ym1D);
    double compute_GCV(const vec & eta, const vec & Y,
                       const vec & ym,  const vec & ym1D);
    void compute_ImAlambda_AYR(const vec & eta, const vec & Y,
                               const vec & ym,  const vec & ym1D);
    
    vec RKHS_theta_sigma(const mat & newX, const vec & ym1D, double sigmaErr);
    void set_dimNull(int dimNull_){
        dimNull = dimNull_;
    }

    virtual mat computeKernel(const mat & dataX1, const mat & dataX2) {
        return Kmat;
    }
    
    virtual mat computeNull(const mat & dataX1) {
        return Smat;
    }
    
    int get_dimNull(){ return dimNull; }
    int get_dimTotal(){ return dimTotal; }
    mat get_Smat(){return Smat; }
    mat get_Kmat(){return Kmat; }
private:
    // gamma is of dimTotal
    int dimNull, dimTotal; // dimension of null space
    
    // n-by-dimNull, n-by-n, (dimNull+n)-by-n
    mat Smat, Kmat, KmatLarge, ImAlambda;
    vec AYR, WYR;
};

typedef std::shared_ptr<thetaRKSH> thetaRKSH_pointer;


#endif
