#include "methodRKHS.hpp"



// ojbective: quadratic loss and penalty
double calModel::objFun(vec gamma){
    double result;
    vec yDiff, theta, thetaNew; // theta to mat
    
    theta = thetaCoreVec[0]->theta_value(gamma);
    thetaNew = thetaLinkVec[0]->theta_trans(theta);
    yDiff = simulator->yModel(thetaNew, dataXReal) - dataY;
    result = 0.5 * sum(yDiff % yDiff);
    result += thetaCoreVec[0]->penalty(gamma) * static_cast<double>(sampleSize);
    
    return result;
}

// gradient of the objective function
vec calModel::gradFun(vec gamma){  // output the gradien
    vec  theta, thetaNew, theta_deriv_multi;
    mat theta_deriv;
    theta = thetaCoreVec[0]->theta_value(gamma);
    theta_deriv = (thetaCoreVec[0]->theta_deriv(gamma));
    thetaNew = thetaLinkVec[0]->theta_trans(theta);
    theta_deriv_multi = thetaLinkVec[0]->theta_trans_deriv(theta);
    
    vec result, yDiff, yPartial;
    yDiff = simulator->yModel(thetaNew, dataXReal) - dataY;
    yPartial = simulator->yModelPartial(thetaNew, dataXReal);
    result =  theta_deriv * (yDiff % yPartial % theta_deriv_multi);
    result += thetaCoreVec[0]->penalty_deriv(gamma) * static_cast<double>(sampleSize);
    return result;
}


// Inverse of Covariance of gamma hat
// Computed from the information matrix
mat calModel::fittedCovMat(vec gamma){
    vec  theta, thetaNew, theta_deriv_multi, theta_deriv_multi2O;
    vec  derivTmp, yHat, yDiff, yPartial, yPartial2O;
    mat theta_deriv;
    
    //deepest level
    theta = thetaCoreVec[0]->theta_value(gamma);
    theta_deriv = thetaCoreVec[0]->theta_deriv(gamma);
    //next level
    theta_deriv_multi = thetaLinkVec[0]->theta_trans_deriv(theta);
    theta_deriv_multi2O = thetaLinkVec[0]->theta_trans_deriv2O(theta); // 2Order Deriv
    //first level
    thetaNew = thetaLinkVec[0]->theta_trans(theta);
    yHat = simulator->yModel(thetaNew, dataXReal);
    yDiff = dataY - yHat;
    yPartial = simulator->yModelPartial(thetaNew, dataXReal);
    yPartial2O = simulator->yModelPartial2O(thetaNew, dataXReal); // 2Order Deriv
    
    derivTmp = yPartial % yPartial % theta_deriv_multi % theta_deriv_multi;
    derivTmp -=  yDiff % yPartial2O % theta_deriv_multi % theta_deriv_multi;
    derivTmp -=  yDiff % yPartial % theta_deriv_multi2O;
    
    mat result;
    result = mat(gamma.n_elem, gamma.n_elem, fill::zeros);
    for(int s = 0; s < sampleSize; s++){
        result += derivTmp(s) * theta_deriv.col(s) * theta_deriv.col(s).t();
    }
    
    vec coef = yDiff % yPartial % theta_deriv_multi;
    result -= thetaCoreVec[0]->theta_deriv_2O(gamma, coef);
    
    result += thetaCoreVec[0]->penalty_deriv2O(gamma) *  static_cast<double>(sampleSize); // 2Order Deriv
    
    // vec actual1O;
    // double residSigmaHat;
    // actual1O = yPartial % theta_deriv_multi;
    // residSigmaHat = thetaCoreVec[0]->compute_residual_sigma(theta, dataY, yHat, actual1O);
    double residSigmaHat = stddev(yDiff);
    result = result / pow(residSigmaHat, 2.0);
    
    
    // actual2O = yPartial % theta_deriv_multi2O +
    //     yPartial2O % theta_deriv_multi % theta_deriv_multi;
    // Rcpp::Rcout << "Sigma 1/2 = " << s1 << " " << s2 << std::endl;
    // result = theta_deriv.t() * result * theta_deriv;
    return result;
}


vec calModel::compute_sigma(const vec & gamma, const mat & newX){
    mat  cchol, gammaCovInv, gammaCov, theta_deriv;
    vec thetaNew, theta_deriv_multi, thetaSigma;
    
    //thetaNew = thetaCoreVec[0]->theta_predict(gamma, newX);
    theta_deriv = thetaCoreVec[0]->theta_deriv_newX(gamma, newX); // at new dataX
    // theta_deriv_multi = thetaLinkVec[0]->theta_trans_deriv(thetaNew);
    // for(int i = 0; i < theta_deriv_multi.n_elem; i++){
    //     theta_deriv.col(i) *= theta_deriv_multi(i);
    // }
    
    // compute the inverse of Sigma (gamma)
    // make it PSD
    vec eigenVal;
    mat eigenVec;
    double invV;
    gammaCovInv = fittedCovMat(gamma);
    eig_sym(eigenVal, eigenVec, gammaCovInv);
    gammaCov = eigenVec.t() * theta_deriv;
    
    for(int i=0; i<eigenVal.n_elem; i++){
        invV = eigenVal(i);
        if(invV < 1e-50){
            invV = 0.0;
        }else{
            invV = 1/sqrt(invV);
        }
        gammaCov.row(i) *= invV;
    }
    
    gammaCov = gammaCov.t() * gammaCov;
    thetaSigma = sqrt(abs(gammaCov.diag()));
    return thetaSigma;
}


vec calModel::adjust_sigma_1O(const vec & sigmaVec, const vec & gradVec){
    vec sigmaAdj(sigmaVec.n_elem);
    sigmaAdj = sigmaVec % abs(gradVec);
    return sigmaAdj;
}


// Second order adjustment of sigma after mapping f(x)
// X: sigma ---> F(x): |f'(x)|sigma
vec calModel::adjust_sigma_2O(const vec & sigmaVec, const vec & gradVec,
                              const vec & hessVec){
    vec sigmaSq;
    vec sigmaAdj(sigmaVec.n_elem);
    sigmaSq = sigmaVec % sigmaVec;
    sigmaAdj =  gradVec % gradVec % sigmaSq +
        0.5 * hessVec % hessVec % sigmaSq % sigmaSq;
    sigmaAdj = sqrt(sigmaAdj);
    return sigmaAdj;
}


// Predict at new value, newX
// With optimized parameter gamma
// Return thetaNew, thetaNewSigma, yHat, yHatSigma, 
mat calModel::predict_y(vec gamma, mat newX, bool getSigma){
    vec theta, yHat, thetaNew, thetaSigma, thetaNewSigma, yHatSigma;
    vec theta_m1, theta_m2, yModel_multiplier;
    vec dirV1, dirV2, yModel_multiplier2;
    mat result;
    
    theta = thetaCoreVec[0]->theta_predict(gamma, newX);
    thetaNew = thetaLinkVec[0]->theta_trans(theta);
    yHat = simulator->yModel(thetaNew, newX);
    if(getSigma){
        if(thetaCoreVec[0]->get_isRKHS()){
            thetaSigma = compute_sigma_version2(gamma, newX);  // numeric stable version for RKHS
        }else{
            // Sigma of H alpha + K gamma
            // In fact, this is the sigma for the parameteric models.
            // The sigma is obtained from the information matrix. 
            thetaSigma = compute_sigma(gamma, newX); 
        }
        theta_m1 = thetaLinkVec[0]->theta_trans_deriv(theta);
        theta_m2 = thetaLinkVec[0]->theta_trans_deriv2O(theta);
        thetaNewSigma = adjust_sigma_2O(thetaSigma, theta_m1, theta_m2);
        
        yModel_multiplier = simulator->yModelPartial(thetaNew, newX);
        yModel_multiplier2 = simulator->yModelPartial2O(thetaNew, newX);
        
        dirV1 = theta_m1 % yModel_multiplier;
        dirV2 = theta_m1 % theta_m1 % yModel_multiplier2 +
            theta_m2 % yModel_multiplier;
        yHatSigma = adjust_sigma_1O(thetaSigma, dirV1);
    }else{
        yHatSigma = vec(newX.n_rows, fill::zeros);
        thetaNewSigma = vec(newX.n_rows, fill::zeros);
        thetaSigma = vec(newX.n_rows, fill::zeros);
        dirV1 = thetaSigma;
        dirV2 = thetaSigma;
        theta_m1 = yHatSigma;
        theta_m2 = yHatSigma;
    }
    
    result = join_horiz(yHat, yHatSigma);
    result = join_horiz(result, thetaNew);
    result = join_horiz(result, thetaNewSigma);
    // result = join_horiz(result, dirV1);
    // result = join_horiz(result, dirV2);
    // 
    return result;
}


double calModel::compute_GCV(vec gamma){
    vec  theta, thetaNew, theta_deriv_multi;
    vec  actual1O, yHat, yPartial;
    double s2;

    theta = thetaCoreVec[0]->theta_value(gamma);
    theta_deriv_multi = thetaLinkVec[0]->theta_trans_deriv(theta);
    thetaNew = thetaLinkVec[0]->theta_trans(theta);
    yHat = simulator->yModel(thetaNew, dataXReal);
    yPartial = simulator->yModelPartial(thetaNew, dataXReal);
    actual1O = yPartial % theta_deriv_multi;
    s2 = thetaCoreVec[0]->compute_GCV(theta, dataY, yHat, actual1O);
    return s2;
}

vec calModel::compute_residualSigma(vec gamma){
    vec theta, thetaNew, theta_deriv_multi;
    vec actual1O, yHat, yPartial;
    vec result(2);
    
    theta = thetaCoreVec[0]->theta_value(gamma);
    theta_deriv_multi = thetaLinkVec[0]->theta_trans_deriv(theta);
    thetaNew = thetaLinkVec[0]->theta_trans(theta);
    yHat = simulator->yModel(thetaNew, dataXReal);
    yPartial = simulator->yModelPartial(thetaNew, dataXReal);
    actual1O = yPartial % theta_deriv_multi;
    result(0) = stddev(dataY - yHat);
    result(1) = thetaCoreVec[0]->compute_residual_sigma(theta, dataY, yHat, actual1O);
    return result;
}



// copmute with direct Bayesian
vec calModel::compute_sigma_version2(const vec & gamma, const mat &newX){
    vec dirV, theta, thetaNew, theta_deriv_multi, thetaSigma;
    theta = thetaCoreVec[0]->theta_value(gamma);
    theta_deriv_multi = thetaLinkVec[0]->theta_trans_deriv(theta);
    thetaNew = thetaLinkVec[0]->theta_trans(theta);
    dirV = simulator->yModelPartial(thetaNew, dataXReal);
    dirV %= theta_deriv_multi;
    
    vec sigmaAll = compute_residualSigma(gamma);
    thetaSigma = thetaCoreVec[0]->RKHS_theta_sigma(newX, dirV, sigmaAll(1)) ; // Sigma of H alpha + K gamma
    return thetaSigma;

}

