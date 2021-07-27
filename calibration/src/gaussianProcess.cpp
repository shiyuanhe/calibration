#include "gaussianProcess.hpp"

gaussianP::gaussianP(){
    nugget = 0.1;
}
// LOO CV Objective function
double gaussianP::cv_obj(vec beta){
    set_beta(beta);
    
    double cv_negLoglik, mu_i, sigmaSq_i, diff;
    
    cv_negLoglik = 0;
    for(int i = 0; i < nSample; i++){
        mu_i = y(i) - alpha(i) / Kinv(i,i);
        sigmaSq_i = Kinv(i, i);
        diff = y(i) - mu_i;
        cv_negLoglik += 0.5 * log(sigmaSq_i) + 0.5 * diff * diff / sigmaSq_i;
    }
    
    return cv_negLoglik;
}


// LOO CV Objective function
vec gaussianP::cv_grad(vec beta){
    set_beta(beta);
    vec grad(2);
    
    for(int j = 0; j < 2; j++){
        grad(j) = - cv_grad_j(beta, j +1);
    }
    
    return grad;
}

double gaussianP::cv_grad_j(vec beta, int j){
    mat Zj, ZjKinv;
    Zj = kernelObj.kernelVar(X, j);
    Zj = Kinv * Zj;
    vec ZjAlpha;
    ZjAlpha = Zj * alpha;
    ZjKinv = Zj * Kinv;
    
    double result, cValue;
    
    result = 0.0;
    for(int i = 0; i< nSample; i++){
        cValue = alpha(i) * ZjAlpha(i);
        cValue -= 0.5 * (1.0 + alpha(i) *  alpha(i) / Kinv(i,i)) * ZjKinv(i,i);
        cValue /= Kinv(i,i);
        result += cValue;    
    }
    
    return result;
}


//yHat value
vec gaussianP::predict_y(mat Xnew){
    vec yHat;
    mat KBasis;
    
    KBasis = kernelObj.kernelCov(Xnew.t(), X);
    yHat = KBasis * alpha;
    
    return yHat;
}

// each row of Xnew is a sample, we need to transpose
mat gaussianP::predict_yCov(mat Xnew){
    mat covMat;
    mat KBasis, matCNew;
    
    KBasis = kernelObj.kernelCov(Xnew.t(), X);
    matCNew = kernelObj.kernelVar(Xnew.t(), 0); //kernelObj.kernelCov(Xnew.t(), X);
    covMat = matCNew - KBasis * solve(K,KBasis).t();
    
    return covMat;
}


vec gaussianP::emul_partial(mat Xnew, int c1){
    vec yHat;
    mat KBasis;
    KBasis = kernelObj.kernelCov_partial_cj(Xnew.t(), X, c1);
    // Rcpp::Rcout << size(KBasis) << std::endl;
    // Rcpp::Rcout << size(alpha) << std::endl;
    yHat = KBasis * alpha;
    // Rcpp::Rcout << size(yHat) << std::endl;
    return yHat;
}

// Used for univariate uncertainty quantification
vec gaussianP::emul_partial2O(mat Xnew, int c1){
    vec yHat;
    mat KBasis;
    KBasis = kernelObj.kernelCov_partial2O_cj(Xnew.t(), X, c1);
    yHat = KBasis * alpha;
    return yHat;
}


// Used for multivariate uncertainty quantification
vec gaussianP::emul_partial2OCross(mat Xnew, int c1, int c2){
    vec yHat;
    mat KBasis;
    if(c1==c2)
        KBasis = kernelObj.kernelCov_partial2O_cj(Xnew.t(), X, c1);
    else
        KBasis = kernelObj.kernelCov_partial2OCross(Xnew.t(), X, c1, c2);
    yHat = KBasis * alpha;
    return yHat;
}

void gaussianP::set_beta(vec beta){
    kernelBeta = beta;
    kernelObj.set_beta(beta);
    K = kernelObj.kernelVar(X, 0);
    K.diag() += nugget;
    Kinv = inv(K);
    alpha = solve(K, y);
}

