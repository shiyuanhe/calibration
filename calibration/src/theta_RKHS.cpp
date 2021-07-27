#include "theta_RKHS.hpp"

void thetaRKSH::set_data(const mat & dataX_,  vec lowerB_, vec upperB_){
    thetaBase::set_data(dataX_, lowerB_, upperB_);
    dimTotal = dimNull + sampleSize;
    Kmat = computeKernel(dataXScaled, dataXScaled);
    Smat = computeNull(dataXScaled);
    if(dimNull > 0){
        KmatLarge = join_horiz(Smat, Kmat);
        KmatLarge = KmatLarge.t();
    }else{
        KmatLarge =  Kmat;
    }
    
}


vec thetaRKSH::theta_value(const vec & gamma){
    vec result;
    
    if(dimNull > 0){
        result = Smat * gamma(span(0, dimNull - 1)) + 
            Kmat * gamma(span(dimNull, dimNull + sampleSize - 1));
    }else{
        result =  Kmat * gamma;
    }
    return result;
}

vec thetaRKSH::theta_predict(const vec & gamma, const mat & newX){
    mat newXScaled = scaleData(newX);
    mat KCov = computeKernel(newXScaled, dataXScaled);
    vec result;
    
    if(dimNull > 0){
        mat SCov = computeNull(newXScaled);
        result = SCov * gamma(span(0, dimNull - 1)) + 
            KCov * gamma(span(dimNull, dimNull + sampleSize - 1));
    }else{
        result =  KCov * gamma;
    }
    
    return result;
}

// Compute residual by
// Y_w^T(I-A(lambda))^2 Y_w / tr(I - A(lambda))
double thetaRKSH::compute_residual_sigma(const vec & eta,
                                         const vec & Y,
                                         const vec & ym, 
                                         const vec & ym1D){
    double tmp1, tmp2, result;
    compute_ImAlambda_AYR(eta, Y, ym, ym1D);
    
    //GCV version (3.26)
    tmp1 = sum(ImAlambda.diag());
    tmp2 = sum(AYR % AYR); // (sqrt(loss2D) % YR % AYR)
    result = tmp2/tmp1;
    
    // REML version (3.30) of Gu Cong
    // tmp2 = sum(WYR % AYR);
    // tmp1 = sampleSize - dimNull;
    // result = tmp2/tmp1;
    return sqrt(result);
}



double thetaRKSH::compute_GCV(const vec & eta,
                              const vec & Y,
                              const vec & ym, 
                              const vec & ym1D){
    double tmp1, tmp2, result;
    
    //GCV 3.23
    compute_ImAlambda_AYR(eta, Y, ym, ym1D);
    tmp1 = sum(ImAlambda.diag()) / static_cast<double>(sampleSize);
    tmp2 = sum(AYR % AYR) / static_cast<double>(sampleSize);
    result = tmp2/(tmp1 * tmp1);
    
    //Gu Cong, (3.29)
    // vec eigVal = eig_sym(ImAlambda);
    // for(int i = 0; i < eigVal.n_elem; i++){
    //     if(eigVal(i) > 1.0e-10){
    //         tmp1 *= eigVal(i);
    //     }
    // }
    // tmp2 = sum(WYR % AYR) / static_cast<double>(sampleSize); 
    // result = tmp2 / pow(tmp1, 1.0/(sampleSize - dimNull));
    return result;
}


// Compute I - A(\lambda) in the paper, stored in . 
// AYR^2/SampleSize in the denominator of GCV 
void thetaRKSH::compute_ImAlambda_AYR(const vec & eta, const vec & Y,
                                      const vec & ym,  const vec & ym1D){
    vec  resid, wiSqrt;
    wiSqrt = abs(ym1D);
    resid = Y - ym;
    
    // Y_w weighted residual
    WYR = eta % wiSqrt - sign(ym1D) % resid;
    
    mat Qw, QRQ, QRR, F2, Qchol;
    Qw = diagmat(wiSqrt) * Kmat * diagmat(wiSqrt);
    if(dimNull > 0){
        mat Sw = diagmat(wiSqrt) * Smat;
        qr(QRQ, QRR, Sw);
        F2 = QRQ.cols(dimNull, sampleSize - 1);
    }else{
        F2 = eye(sampleSize, sampleSize);
    }
    
    ImAlambda = F2.t() * Qw * F2;
    ImAlambda.diag() +=  2.0 * sampleSize * lambda;
    
    chol(Qchol, ImAlambda);
    F2 = solve(trimatl(Qchol.t()), F2.t());
    
    ImAlambda = (2.0 * sampleSize * lambda) * F2.t()  * F2; // why times 2
    
    // Residual = Y_w - A(\lambda) Y_w for the weighted response.
    AYR = ImAlambda * WYR;
    
}



vec thetaRKSH::RKHS_theta_sigma(const mat & newX,
                                const vec & ym1D,
                                double sigmaErr){
    mat Mmat, Mmat_R, Sw;
    vec wiSqrt = abs(ym1D), phi, xi, tVec;
    int i;
    double tmp1, tmp2, bDbl, rho;
    bDbl = sigmaErr * sigmaErr / 
        (2.0 * static_cast<double>(sampleSize) * lambda);
    
    rho = 5e4;
    Sw = Smat; //diagmat(wiSqrt)
    Mmat =  Kmat;
    Mmat.diag() += (2.0 * static_cast<double>(sampleSize) * lambda) / (1e-20 + wiSqrt % wiSqrt);
    if(dimNull > 0)
        Mmat += rho * Sw *Sw.t();
    chol(Mmat_R, Mmat);
    
    mat newXScaled = scaleData(newX);
    mat KmatNew = computeKernel(newXScaled, newXScaled);
    mat SmatNew = computeNull(newXScaled);
    mat KCov = computeKernel(newXScaled, dataXScaled);
    vec result(newX.n_rows, fill::zeros);
    for(i = 0; i < newX.n_rows; i++){
        xi = KCov.row(i).t();
        tmp1 = KmatNew(i,i);
        tVec = xi;
        if(dimNull > 0){
            phi = SmatNew.row(i).t();
            tmp1 += rho * sum(phi % phi);
            tVec += rho * Sw * phi;
        }
        tVec = solve(trimatl(Mmat_R.t()), tVec);
        tmp2 = sum(tVec % tVec);
        result(i) = (tmp1 - tmp2) * bDbl;
    }
    
    return sqrt(result);
}

