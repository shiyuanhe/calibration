#include "methodBayes.hpp"

mcmcCal::mcmcCal(int problemIndex, int thetaIndex, int linkIndex, int numComp) : 
    methodBase(problemIndex, thetaIndex, linkIndex, numComp)
{
    outSampleRequired = false;
    seqN = 1;
    hasPrior = false;
    
    if (numComp > 1)
    {
        // select objective function for the multivariate case
        switch (problemIndex)
        {
        case 1:
            simulator = std::make_shared<simuobj_multi1>();
            break;
        case 2:
            simulator = std::make_shared<simuobj_multi2>();
            break;
        default:
            simulator = std::make_shared<simuobj_GP>(); break;
        }
    }
}


void mcmcCal::setPrior(double ay_, double by_,
              double atheta_, double btheta_,
              double mu_theta_)
{
    // set prior
    ay = ay_;
    by = by_;
    atheta = atheta_;
    btheta = btheta_;
    mu_theta = mu_theta_;
    hasPrior = true;
}

// Make prediction outside training set.
// Usually at a dense grid of X
void mcmcCal::setPredition(mat dataXSeq_)
{
    outSampleRequired = true;
    dataXSeq = dataXSeq_;
    seqN = dataXSeq.n_rows;
}

// Main Function
List mcmcCal::execute_mcmc(int nBurnin, int nSample)
{
    int mcmcI;
    initRecordSample(nSample);

    for (mcmcI = 1; mcmcI <= nBurnin; mcmcI++)
    {
        oneRoundSample(false);
        adjustStepSize(mcmcI);
    }

    for (mcmcI = 0; mcmcI < nSample; mcmcI++)
    {
        oneRoundSample(true);
        recordSample(mcmcI);
    }

    return List::create(Named("thetaS") = thetaSamples,
                        Named("yHat") = yHatAll,
                        Named("nu") = nuSamples,
                        Named("lambdaY") = lySamples,
                        Named("lambdaT") = ltSamples,
                        Named("thetaSeqAll") = thetaSeqAll,
                        Named("ySeqAll") = ySeqAll,
                        Named("Count1") = Count1,
                        Named("Count2") = Count2);
}

void mcmcCal::oneRoundSample(bool recordStage)
{
    int compI, iterSub;

    for (compI = 0; compI < numComp; compI++)
    {
        sample_nu(compI); // this is where the eigen decom is created
        // propose thetaVec 5 times more often
        for (iterSub = 0; iterSub < 5; iterSub++)
            sample_thetaFun(compI);
        sample_lambda_theta(compI);
        if (outSampleRequired && recordStage)
            compute_thetaSeq(compI);
    }
    yHat = compute_yHat(thetaVec);
    sample_lambda_y();
}

void mcmcCal::adjustStepSize(int mcmcI)
{
    // Adjust step size by acceptance rate, for better mixing
    // Like the trust region method
    if (mcmcI % 200 == 0)
    {
        if (Count1 > 130)
            proposeC1 *= 1.5;
        if (Count1 < 100)
            proposeC1 *= 0.5;
        if (Count2 / numComp > 600)
            proposeC2 *= 1.5;
        if (Count2 / numComp < 400)
            proposeC2 *= 0.5;
        Count1 = 0;
        Count2 = 0;
    }
}

// This init function shold be rewritten
void mcmcCal::init_mcmc(double proposeC1_, double proposeC2_,
                        double nu_, double lambda_y_,
                        double lambda_theta_)
{
    if(sampleSize == 0) throw(std::runtime_error("Please set sample first!"));
    if(!hasPrior) throw(std::runtime_error("Pior not specified yet!"));
    
    // stepsize for nu
    proposeC1 = proposeC1_;
    // stepsize for theta
    proposeC2 = proposeC2_;
    
    // For KernelFun, the input should be stored by column order
    vec scaling_ = ones(dataXReal.n_cols);
    kernelFun.set_scaling(scaling_);
    dataXRealTrans = dataXReal.t();
    dataXSeqTrans = dataXSeq.t();
    
    nu = ones(numComp) * nu_;
    thetaVec = ones(sampleSize, numComp) * mu_theta; //mu_theta * ones(sampleSize);
    thetaVecNew = mat(sampleSize, numComp);
    lambda_y = lambda_y_;
    lambda_theta = arma::ones(numComp) * lambda_theta_;
    
    Lambda.resize(numComp);
    U.resize(numComp);
    // Random perturb thetaVecInit by a large amount.
    // The input sample is initialized at the prior mean.
    // Large perturbation is performed for random init.
    for (int compI = 0; compI < numComp; compI++)
    { 
        mat Rnu = compute_kernel(nu(compI));
        vec testValue;
        mat testVec;
        eig_sym(testValue, testVec, Rnu);
        Lambda.at(compI) = testValue;
        U.at(compI) = testVec;
        proposeC2 = 0.5;
        thetaVec = propose_thetaFun(compI);
        thetaVecNew.col(compI) = thetaLinkVec[compI]->theta_trans(thetaVec.col(compI));
    }
    yHat = compute_yHat(thetaVec); // accept any way
    
    // set back propose step size
    proposeC1 = proposeC1_;
    proposeC2 = proposeC2_;
    // Warm up by first round of MCMC
    oneRoundSample(false);

    // Init acceptance rate
    Count2 = 0;
    Count1 = 0;
}

// Init proper container size to record MCMC samples
void mcmcCal::initRecordSample(int nSample)
{
    yHatAll = mat(sampleSize, nSample);
    ySeqAll = mat(seqN, nSample);
    nuSamples = mat(numComp, nSample);
    ltSamples = mat(numComp, nSample);
    lySamples = vec(nSample);
    thetaSamples = cube(sampleSize, numComp, nSample);
    thetaSeqAll = cube(seqN, numComp, nSample);
    thetaSeq  = mat(seqN, numComp);
}

// record current MCMC samples
void mcmcCal::recordSample(int mcmcI)
{
    yHatAll.col(mcmcI) = yHat;
    nuSamples.col(mcmcI) = nu;
    lySamples(mcmcI) = lambda_y;
    ltSamples.col(mcmcI) = lambda_theta;
    thetaSamples.slice(mcmcI) = thetaVecNew;

    // sample at a dense new location
    if (outSampleRequired)
    {
        compute_ySeq();
        thetaSeqAll.slice(mcmcI) = thetaSeq;
        ySeqAll.col(mcmcI) = ySeq;
    }
}

// Sample new lambda_y
// It depends on prior ay, by and yHat
void mcmcCal::sample_lambda_y()
{
    //Rcpp::Rcout << "sample_lambda_y" << std::endl;

    double gamma_a, gamma_b;
    gamma_a = ay + static_cast<double>(sampleSize) / 2.0;
    gamma_b = 1.0 / by - compute_expS1(yHat); // scale parameter
    lambda_y = rgamma(1, gamma_a, 1.0 / gamma_b)[0];
}


// Required U/Lambda decompostive of R_nu
void mcmcCal::sample_lambda_theta(int compI)
{
    //Rcpp::Rcout << "sample_lambda_theta" << std::endl;

    double gamma_a, gamma_b;
    gamma_a = atheta + static_cast<double>(sampleSize) / 2.0;
    gamma_b = 1.0 / btheta - 
        compute_expS2(thetaVec.col(compI), U.at(compI), Lambda.at(compI)); // scale parameter
    lambda_theta(compI) = rgamma(1, gamma_a, 1.0 / gamma_b)[0];
}

// Update nu for component, compI
// Required: previous U/Lambda, nu, lambda_theta
void mcmcCal::sample_nu(int compI)
{
    //Rcpp::Rcout << "sample_nu" << std::endl;
    double nu_propose, likRatio, likS1, likS2;
    mat R_nu_propose, U_propose;
    vec Lambda_propose;

    // propose new nu for compI, and compute the correlation matrix.
    // The eigen decomposition is performed for convenience
    // of other related functions.
    nu_propose = nu(compI) + rnorm(1, 0, proposeC1)[0];
    
    R_nu_propose = compute_kernel(nu_propose);
    eig_sym(Lambda_propose, U_propose, R_nu_propose);
    
    //Lambda_propose += 1e-2;
    // Gamble forward
    likS1 = computeLikS_nu(nu(compI), U.at(compI), Lambda.at(compI), compI);
    likS2 = computeLikS_nu(nu_propose, U_propose, Lambda_propose, compI);
    likRatio = exp(likS2 - likS1);
    if (runif(1, 0, 1)[0] < likRatio)
    {
        Count1++;
        nu(compI) = nu_propose;
        U.at(compI) = U_propose;
        Lambda.at(compI) = Lambda_propose;
    }
}

// compute the correlation matrix
// rho = exp( -exp(nu) )
mat mcmcCal::compute_kernel(double nu_)
{
    mat R_nu_propose;
    vec kernel_beta = zeros(2);
    kernel_beta(1) = -nu_; //exp(nu_) - 2;
    kernelFun.set_beta(kernel_beta);
    R_nu_propose = kernelFun.kernelVar(dataXRealTrans, 0);
    R_nu_propose.diag() += 1e-8;
    return R_nu_propose;
}

// Compute the loglik for the proposed nu_, and
// the corresponding eigen decompostion of R_nu
double mcmcCal::computeLikS_nu(double nu_, const mat &U_,
                               const vec Lambda_, const int compI)
{
    double tmp, result;
    result = lambda_theta(compI) * 
        compute_expS2(thetaVec.col(compI), U_, Lambda_);
    result += nu_ - exp(nu_);
    result -= 0.5 * sum(log(Lambda_));
    tmp = 1 - exp(-exp(nu_)); // ?? nu_ or nu
    result += (btheta - 1) * log(tmp);
    return result;
}

// Update for the compI-th column.
// U, Lambda should be computed before this.
void mcmcCal::sample_thetaFun(int compI)
{
    //Rcpp::Rcout << "sample_thetaVec" << std::endl;

    double likS1, likS2, likRatio;
    mat thetaVec_propose;
    vec yHat_propose;

    thetaVec_propose = propose_thetaFun(compI);
    yHat_propose = compute_yHat(thetaVec_propose);
    
    // evaluate the likelihood ratio for the I-th component
    likS1 = lambda_y * compute_expS1(yHat) +
            lambda_theta(compI) * 
            compute_expS2(thetaVec.col(compI), U.at(compI), Lambda.at(compI));
    
    // Is likS2 computed correctly??
    likS2 = lambda_y * compute_expS1(yHat_propose) +
            lambda_theta(compI) *
                compute_expS2(thetaVec_propose.col(compI), U.at(compI), Lambda.at(compI));

    likRatio = exp(likS2 - likS1);
    if (runif(1, 0, 1)[0] < likRatio)
    {
        Count2++;
        yHat = yHat_propose;
        thetaVec.col(compI) = thetaVec_propose.col(compI);
        // thetaVecNew is not updated when computing yHat
        thetaVecNew.col(compI) = thetaLinkVec[compI]->theta_trans(thetaVec.col(compI));
        
    }
}

// Update for the compI-th column.
// U, Lambda should be computed before this.
// Return proposed theta vector, only for the specified component.
mat mcmcCal::propose_thetaFun(int compI)
{
    mat thetaVec_propose;
    vec addZ;
    thetaVec_propose = thetaVec;
    addZ = vec(rnorm(sampleSize));
    addZ %= (proposeC2 * sqrt(Lambda.at(compI)));
    thetaVec_propose.col(compI) += U.at(compI) * addZ;
    return thetaVec_propose;
}


// Compute model output, yHat, from  untransformed theta.
vec mcmcCal::compute_yHat(const mat &theta_)
{
    mat thetaTrans(sampleSize, numComp);
    vec yHat_;
    for (int i = 0; i < numComp; i++)
        thetaTrans.col(i) = thetaLinkVec[i]->theta_trans(theta_.col(i));
    yHat_ = simulator->yModel(thetaTrans, dataXReal);
    return yHat_;
}

// Sampling at new position
void mcmcCal::compute_ySeq()
{
    // transform and compute yHat
    for (int i = 0; i < numComp; i++)
        thetaSeq.col(i) = thetaLinkVec[i]->theta_trans(thetaSeq.col(i));

    ySeq = simulator->yModel(thetaSeq, dataXSeq);
}

void mcmcCal::compute_thetaSeq(int compI)
{
    mat Rcov, RVar, RMinus, RChol, tmp;
    vec kernel_beta = zeros(2);

    kernel_beta(1) = -nu(compI); //exp(nu_) - 2;
    kernelFun.set_beta(kernel_beta);

    // conditional mean, thetaSeq given thetaVec
    Rcov = kernelFun.kernelCov(dataXSeqTrans, dataXRealTrans);
    tmp = (U.at(compI).t() * (thetaVec.col(compI) - mu_theta)) / Lambda.at(compI);
    thetaSeq.col(compI) = Rcov * U.at(compI) * tmp;

    // conditional variance
    RVar = kernelFun.kernelVar(dataXSeqTrans, 0);
    RMinus = Rcov * U.at(compI) * diagmat(1 / sqrt(Lambda.at(compI)));
    RVar -= RMinus * RMinus.t();
    RVar.diag() += 0.008;
    RChol = chol(RVar);

    // add perturbation to the conditional mean
    thetaSeq.col(compI) += RChol.t() * vec(rnorm(seqN));
}
