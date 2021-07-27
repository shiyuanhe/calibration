#include "method_RKHS_multivariate.hpp"


calModelMulti::calModelMulti(int c1, int c2, int c3, int numComp_): 
    methodBase(c1, c2, c3, numComp_)
{
    numComp = numComp_;
    sampleSize = 0;
    thetaCoreVec.resize(numComp);
    thetaLinkVec.resize(numComp);
    
    // objective function
    switch (c1)
    {
    case 1:
        simulator = std::make_shared<simuobj_multi1>();
        break;
    case 2:
        simulator = std::make_shared<simuobj_multi2>();
        break;
    default:
        simulator = std::make_shared<simuobj_GP>(); 
    break;
    }
    
}

// ojbective: quadratic loss and penalty
double calModelMulti::objFun(mat gamma)
{
    double result;
    vec yDiff;
    //gamma.reshape(gamma.size()/numComp, numComp);
    
    yDiff = computeYHat(gamma) - dataY;
    result = 0.5 * sum(yDiff % yDiff);
    // add penalty
    double penaltyValue;
    for (size_t i = 0; i < numComp; i++)
    {
        penaltyValue = thetaCoreVec[i]->penalty(gamma.col(i)) * 
            static_cast<double>(sampleSize);
        result += penaltyValue;
    }
    
    return result;
}

// gradient of the objective function
mat calModelMulti::gradFun(mat gamma)
{
    mat theta_deriv; 
    mat result(arma::size(gamma)), yPartial;
    vec yDiff;
    
    yDiff = computeYHat(gamma) - dataY;
    yPartial = computeYHat_1O(gamma);
    //gamma.reshape(gamma.size()/numComp, numComp);
    
    // Theta: parameter number X # of observations
    for (size_t i = 0; i < numComp; i++)
    {
        theta_deriv = thetaCoreVec[i]->theta_deriv(gamma.col(i));
        result.col(i) = theta_deriv * (yDiff % yPartial.col(i));
        result.col(i) += thetaCoreVec[i]->penalty_deriv(gamma.col(i)) * 
            static_cast<double>(sampleSize);
    }
    
    return result;
}

// compute the transformed theta from parameters gamma
mat calModelMulti::computeTheta(const mat  & gamma){
    mat theta(sampleSize, numComp), thetaNew(sampleSize, numComp); // theta to mat
    size_t i;
    
    // Compute Theta by direct RKHS
    // Possibly transform it to proper range
    for (i = 0; i < numComp; i++)
    {
        theta.col(i) = thetaCoreVec[i]->theta_value(gamma.col(i));
        thetaNew.col(i) = thetaLinkVec[i]->theta_trans(theta.col(i));
    }
    return thetaNew;
}


// Compute the response $y^p\circ l$ from the raw gamma
// The effect of the link is included.
vec calModelMulti::computeYHat(const mat & gamma){
    mat thetaNew = computeTheta(gamma);
    vec yhat = simulator->yModel(thetaNew, dataXReal);
    return yhat;
}


// Compute the first order derivative of the model
// $y^p\circ l$, the effect of the link is included.
// Returned Matrix: # of rows is the sample size, # of cols is the comp number. 
mat calModelMulti::computeYHat_1O(const mat & gamma){
    mat theta(sampleSize, numComp), thetaNew(sampleSize, numComp);
    // The multiplicative factor due to tranformation
    mat theta_deriv_multi(sampleSize, numComp);        
    
    size_t i;
    for (i = 0; i < numComp; i++){
        theta.col(i) = thetaCoreVec[i]->theta_value(gamma.col(i));
        thetaNew.col(i) = thetaLinkVec[i]->theta_trans(theta.col(i));
        theta_deriv_multi.col(i) = thetaLinkVec[i]->theta_trans_deriv(theta.col(i));
    }
    
    mat  yPartial;
    yPartial = simulator->yModelPartial(thetaNew, dataXReal);
    
    return (yPartial % theta_deriv_multi);
}

// Predict at new value, newX
// With optimized parameter gamma
mat calModelMulti::predict_y(mat gamma, mat newX, bool getSigma)
{
    mat result(newX.n_rows, numComp + 1);
    
    // The multiplicative factor due to tranformation
    mat theta(newX.n_rows, numComp);
    mat thetaNew(newX.n_rows, numComp); // theta to mat
    size_t i;
    
    // Compute Theta by direct RKHS
    // Possibly transform it to proper range
    for (i = 0; i < numComp; i++)
    {
        theta.col(i) = thetaCoreVec[i]->theta_predict(gamma.col(i),  newX);
        thetaNew.col(i) = thetaLinkVec[i]->theta_trans(theta.col(i));
    }
    
    
    result.col(0) = simulator->yModel(thetaNew, newX);
    result.cols(1, numComp) = thetaNew;
    return result;
}


// Compute the first order derivative at new value, newX
// with optimized parameter gamma
mat calModelMulti::predict_yDeriv(mat gamma, mat newX)
{
    // The multiplicative factor due to tranformation
    mat theta_deriv_multi(newX.n_rows, numComp);
    mat theta(newX.n_rows, numComp);
    mat thetaNew(newX.n_rows, numComp); // theta to mat
    size_t i;
    
    // Compute Theta by direct RKHS
    // Possibly transform it to proper range
    for (i = 0; i < numComp; i++)
    {
        theta.col(i) = thetaCoreVec[i]->theta_predict(gamma.col(i),  newX);
        thetaNew.col(i) = thetaLinkVec[i]->theta_trans(theta.col(i));
        theta_deriv_multi.col(i) = thetaLinkVec[i]->theta_trans_deriv(theta.col(i));
    }
    
    
    mat  yPartial;
    yPartial = simulator->yModelPartial(thetaNew, newX);
    return (yPartial % theta_deriv_multi);
}

// Compute the generalized CV.
double calModelMulti::compute_GCV(mat gamma, double lambda){
    vec tmp = compute_IminusA(gamma, lambda);
    return tmp(0) / (tmp(1) * tmp(1));
}


// Compute the residual sigma sqaured
double calModelMulti::compute_sigmaSq(mat gamma, double lambda){
    vec tmp = compute_IminusA(gamma, lambda);
    return tmp(0) / tmp(1);
}


vec calModelMulti::compute_IminusA(mat gamma, double lambda){
    
    mat thetaNew, yHatDeriv;
    vec yHat, yHat_weighted;
    
    thetaNew = computeTheta(gamma);
    yHat = computeYHat(gamma);
    yHatDeriv = computeYHat_1O(gamma);
    
    // Compute the new response as Y - Yhat + \sum w_i theta_i
    // w_i is the first order derivative of the model w.r.t theta_i
    size_t i,j;
    yHat_weighted = dataY - yHat;
    for(i = 0; i < numComp; i++){
        yHat_weighted += yHatDeriv.col(i) % thetaNew.col(i);
    }
    
    // Get the null and native dimension of each RKHS
    vec dimNull(numComp), dimRKHS(numComp);
    thetaRKSH_pointer castedPointer;
    for(i = 0; i < numComp; i++){
        castedPointer = std::dynamic_pointer_cast<thetaRKSH>(thetaCoreVec[i]);
        dimNull(i) = castedPointer -> get_dimNull();
        dimRKHS(i) = castedPointer -> get_dimTotal() - dimNull(i);
    }
    
    // Compute the stacked response yTotal, null matrix Stotal,
    // and kernel matrix Ktotal.
    mat Stotal(numComp * sampleSize, sum(dimNull));
    mat Ktotal(numComp * sampleSize, sum(dimRKHS));
    vec yTotal(numComp * sampleSize);
    
    mat tmpS, tmpK;
    int SColstartI, KColstartI;
    arma::span rowSpan, colSpan1, colSpan2;
    
    SColstartI = 0;
    KColstartI = 0;
    
    for(i = 0; i < numComp; i++){
        castedPointer = std::dynamic_pointer_cast<thetaRKSH>(thetaCoreVec[i]);
        
        // Store them as the weighted version
        tmpS = castedPointer -> get_Smat();
        tmpK = castedPointer -> get_Kmat();
        
        tmpS.each_col() %= yHatDeriv.col(i);
        tmpK.each_row() %= yHatDeriv.col(i).t();
        tmpK.each_col() %= yHatDeriv.col(i);
        
        // Put them into appropriate position
        // Duplicate each component, put it in the rows
        colSpan1 = span(SColstartI, SColstartI + dimNull(i) - 1);
        colSpan2 = span(KColstartI, KColstartI + dimRKHS(i) - 1);
        
        for(j = 0; j < numComp; j++){
            rowSpan = span(j * sampleSize, (j+1) * sampleSize - 1);
            if(i == 0) yTotal(rowSpan) = yHat_weighted;
            // update null matrix
            Stotal(rowSpan, colSpan1) = tmpS;
            // update kernel matrix
            Ktotal(rowSpan, colSpan2) = tmpK;
        }
        SColstartI += dimNull(i);
        KColstartI += dimRKHS(i);
    }
    
    
    // Apply QR decomposition to Stotal
    mat QR_Q, QR_R, F2;
    qr(QR_Q, QR_R, Stotal);
    F2 = QR_Q.cols(sum(dimNull), numComp * sampleSize - 1);
    
    mat FPhiF_I;
    FPhiF_I = F2.t() * Ktotal * F2;
    FPhiF_I.diag() += 2.0 * sampleSize * lambda;
    
    mat IminusA;
    IminusA = solve(FPhiF_I, F2.t());
    IminusA = F2 * IminusA;
    IminusA *= 2.0 * sampleSize * lambda;
    
    double numerGCV, denomGCV;
    yTotal = IminusA * yTotal;
    denomGCV = sum(IminusA.diag());
    numerGCV = sum(yTotal % yTotal);
    numerGCV /= static_cast<double>(sampleSize);
    denomGCV /= static_cast<double>(sampleSize);
    numerGCV /= static_cast<double>(numComp);
    denomGCV /= static_cast<double>(numComp);
    
    vec res(2);
    res(0) = numerGCV;
    res(1) = denomGCV;
    
    return res;
}



// Compute the the prediction error
mat calModelMulti::sigmaSq_y(mat gamma, mat newX, double lambda){
    mat result;
    if(thetaCoreVec[0]->get_isRKHS()){
        result = RKHS_sigmaY(gamma, newX, lambda);
    }else{
        result = Para_sigmaY(gamma, newX, lambda);
    } 
    // we actually obtained variance in the previous step.
    result = sqrt(result);
    return result;
}



// Compute the generalized CV.
mat calModelMulti::sigmaSq_y_simu(mat newX, mat yNew, mat ySigma, double err_epsilon){
    int monteN = 2000;
    mat result(newX.n_rows, monteN);
    mat thetaMonte(monteN, numComp);
    mat NewXDup;
    if(thetaCoreVec[0]->get_isRKHS()){
        
        int i;
        for(i = 0; i < newX.n_rows; i++){
            
            //Rcpp::Rcout << ySigma(i,0) / abs(yNew(i,0))  << endl;
            if(ySigma(i,0) / abs(yNew(i,0))  > err_epsilon)
                continue;
            NewXDup = repmat(newX.row(i), monteN, 1);
            thetaMonte = randn(monteN, numComp);
            thetaMonte.each_row() %= ySigma(i,span(1,numComp))*0.1;
            thetaMonte.each_row() += yNew(i,span(1,numComp));
            result.row(i) = (simulator->yModel(thetaMonte, NewXDup)).t();
        }
        
    } 
    return result;
}

// return: ysigma, theta1_sigma, theta2_sigma
mat calModelMulti::RKHS_sigmaY(mat gamma, mat newX, double lambda){
    mat  yHatDeriv, yHatDerivNew;
    mat result(newX.n_rows,numComp+1);
    yHatDeriv = computeYHat_1O(gamma);
    yHatDerivNew = predict_yDeriv(gamma, newX);
    
    // Get the null and native dimension of each RKHS
    mat newXScaled, oldXScaled;
    thetaRKSH_pointer castedPointer;
    castedPointer = std::dynamic_pointer_cast<thetaRKSH>(thetaCoreVec[0]);
    newXScaled = castedPointer->scaleData(newX);
    oldXScaled = castedPointer->scaleData(dataXReal);
    // covariance between yNew and yOld
    mat covSum(newX.n_rows, sampleSize, fill::zeros);
    // covariance for yOld
    mat KwSum(sampleSize, sampleSize, fill::zeros);
    // covariance for yNew
    mat KwSumNew(newX.n_rows, newX.n_rows, fill::zeros);
    
    mat tmpS, tmpK, tmpSNew, tmpKNew, tmpKCov;
    double rho = 5e4;
    size_t i;
    for(i = 0; i < numComp; i++){
        
        castedPointer = std::dynamic_pointer_cast<thetaRKSH>(thetaCoreVec[i]);
        
        // Kernel for old data
        tmpS = castedPointer -> get_Smat();
        tmpK = castedPointer -> get_Kmat();
        tmpK = rho * tmpS * tmpS.t() + tmpK;
        tmpK.each_row() %= yHatDeriv.col(i).t();
        tmpK.each_col() %= yHatDeriv.col(i);
        KwSum += tmpK;  // Lambda_{22} in the paper
        
        // Kernel for new data
        tmpSNew = castedPointer->computeNull(newXScaled);
        //tmpSNew.each_col() %= yHatDerivNew.col(i);
        tmpKNew = castedPointer->computeKernel(newXScaled, newXScaled);
        tmpKNew += rho * tmpSNew * tmpSNew.t();
        tmpKNew.each_col() %= yHatDerivNew.col(i);
        tmpKNew.each_row() %= yHatDerivNew.col(i).t();
        KwSumNew += tmpKNew;
        
        
        // kernel between new and old
        tmpKCov = castedPointer->computeKernel(newXScaled, oldXScaled);
        tmpKCov += rho * tmpSNew * tmpS.t();
        tmpKCov.each_col() %= yHatDerivNew.col(i);
        tmpKCov.each_row() %= yHatDeriv.col(i).t();
        covSum += tmpKCov;
        
    }
    
    
    double sigmaSqRes;
    sigmaSqRes = compute_sigmaSq(gamma, lambda);
    KwSum.diag() += sampleSize *  lambda;
    
    mat tmpResult;
    tmpResult = KwSumNew - covSum * KwSum.i() * covSum.t();
    result.col(0) = tmpResult.diag();
    //yHatDeriv = join_vert(yHatDeriv, yHatDerivNew);
    
    
    mat theta_deriv_multi(newX.n_rows, numComp);
    mat theta(newX.n_rows, numComp);
    mat thetaNew(newX.n_rows, numComp); // theta to mat

    // Compute Theta by direct RKHS
    // Possibly transform it to proper range
    for (i = 0; i < numComp; i++)
    {
        theta.col(i) = thetaCoreVec[i]->theta_predict(gamma.col(i),  newX);
        thetaNew.col(i) = thetaLinkVec[i]->theta_trans(theta.col(i));
        theta_deriv_multi.col(i) = thetaLinkVec[i]->theta_trans_deriv(theta.col(i));
    }
    
    // repeat again to compute for each component
    for(i = 0; i < numComp; i++){
        castedPointer = std::dynamic_pointer_cast<thetaRKSH>(thetaCoreVec[i]);
        tmpS = castedPointer -> get_Smat();
        // Kernel for new data
        tmpSNew = castedPointer->computeNull(newXScaled);
        tmpKNew = castedPointer->computeKernel(newXScaled, newXScaled);
        tmpKNew += rho * tmpSNew * tmpSNew.t();
        //tmpKNew.each_col() = theta_deriv_multi.col(i);
        //tmpKNew.each_row() = theta_deriv_multi.col(i).t();
        
        // kernel between new and old
        tmpKCov = castedPointer->computeKernel(newXScaled, oldXScaled);
        tmpKCov += rho * tmpSNew * tmpS.t();
        //tmpKCov.each_col() %= theta_deriv_multi.col(i);
        tmpKCov.each_row() %= yHatDeriv.col(i).t();
        
        tmpResult = tmpKNew - tmpKCov * KwSum.i() * tmpKCov.t();
        result.col(i+1) = tmpResult.diag() % theta_deriv_multi.col(i) % theta_deriv_multi.col(i);
        

    }
    
    
    return result* (sigmaSqRes / (sampleSize *  lambda));
}


vec calModelMulti::Para_sigmaY(mat gamma, mat newX, double lambda){
    mat  inforMat, yHatDerivNew;
    mat result, theta_deriv, tmp;
    yHatDerivNew = predict_yDeriv(gamma, newX);
    size_t dimP, i;
    dimP = gamma.n_rows;
    theta_deriv = zeros<mat>(dimP * numComp, newX.n_rows);
    for(i = 0; i < numComp; i++){
        tmp = thetaCoreVec[i]->theta_deriv_newX(gamma.col(i), newX);
        tmp.each_row() %= yHatDerivNew.col(i).t();
        theta_deriv.rows(i * dimP, (i+1) * dimP - 1) = tmp;
    }

    inforMat = parametric_info(gamma);
    
    vec eigenVal;
    mat eigenVec;
    eig_sym(eigenVal, eigenVec, inforMat);
    
    for(int i=0; i<eigenVal.n_elem; i++){
        if(eigenVal(i) < 1e-10){
            eigenVal(i) = 0.0;
        }else{
            eigenVal(i) = 1.0/sqrt(eigenVal(i));
        }
    }
    eigenVec.each_row() %= eigenVal.t();
    theta_deriv = eigenVec.t() * theta_deriv;
    
    result = theta_deriv.t() * theta_deriv;
    return result.diag();
}

// Parametric Information matrix
mat calModelMulti::parametric_info(const mat& gamma){
    vec  derivTmp, yHat, yDiff, yPartial2O;
    mat  yPartial, theta_derivI, theta_derivJ, tmp, result;
    
    mat theta_deriv_multi(sampleSize, numComp);
    mat theta_deriv_multi2O(sampleSize, numComp);
    mat theta(sampleSize, numComp);
    mat thetaNew(sampleSize, numComp); // theta to mat
    size_t i,j, dimP;
    

    // theta, link transformation, link derivatives
    for (i = 0; i < numComp; i++)
    {
        // deepest level
        theta.col(i) = thetaCoreVec[i]->theta_predict(gamma.col(i),  dataXReal);
        thetaNew.col(i) = thetaLinkVec[i]->theta_trans(theta.col(i));
        // link level
        theta_deriv_multi.col(i) = thetaLinkVec[i]->theta_trans_deriv(theta.col(i));
        theta_deriv_multi2O.col(i) = 
            thetaLinkVec[i]->theta_trans_deriv2O(theta.col(i));
    }

    yHat = simulator->yModel(thetaNew, dataXReal);
    yDiff = yHat - dataY;
    // w.r.t theta_1, theta_2 ....
    yPartial = simulator->yModelPartial(thetaNew, dataXReal);

    dimP = gamma.n_rows;
    result = zeros<mat>(dimP * numComp, dimP * numComp);
    arma::span sp1, sp2;
    for(i = 0; i < numComp; i++){
        for(j = i; j < numComp; j++){
            // Compute for theta_i and theta_j
            yPartial2O = simulator->yModelPartial2OCross(thetaNew, dataXReal, i, j); 
            derivTmp = yPartial.col(i) % yPartial.col(j) % 
                theta_deriv_multi.col(i) % theta_deriv_multi.col(j);
            derivTmp -=  yDiff % yPartial2O % 
                theta_deriv_multi.col(i) % theta_deriv_multi.col(j);
            if(i==j){
                derivTmp -=  yDiff % yPartial.col(i) % theta_deriv_multi2O.col(i);
            }

            // Sample stored column-wise 
            theta_derivI = thetaCoreVec[i]->theta_deriv(gamma.col(i));
            theta_derivJ = thetaCoreVec[j]->theta_deriv(gamma.col(j));
            theta_derivI.each_row() %= derivTmp.t();
            tmp = theta_derivI * theta_derivJ.t();

            if(i==j){
                vec coef = yDiff % yPartial.col(i) % theta_deriv_multi.col(i);
                tmp -= thetaCoreVec[i]->theta_deriv_2O(gamma.col(i), coef);
            }
            
            // put tmp into result
            sp1 = arma::span(i * dimP, (i+1) * dimP - 1);
            sp2 = arma::span(j * dimP, (j+1) * dimP - 1);
            result(sp1, sp2) = tmp;
            result(sp2, sp1) = tmp;
            
        }

    }
    

    double residSigmaHat = stddev(yDiff);
    result = result / pow(residSigmaHat, 2.0);
    
    return result;
}

