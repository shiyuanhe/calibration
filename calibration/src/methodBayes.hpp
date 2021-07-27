#ifndef CLASS_MCMCCAL
#define CLASS_MCMCCAL

#include "methodBase.hpp"
#include "simuObj_multivariate.hpp"

class mcmcCal : public methodBase
{
  public:
    mcmcCal(int c1, int c2, int c3, int c4);

    void setPrior(double ay_, double by_,
                  double atheta_, double btheta_,
                  double mu_theta_);

    // Make prediction outside training set.
    // Usually at a dense grid of X
    void setPredition(mat dataXSeq_);

    List execute_mcmc(int nBurnin, int nSample);
    void init_mcmc(double, double, double, double, double);

  private:
    squaredExp kernelFun;

    // MCMC containers
    cube thetaSamples, thetaSeqAll;
    mat yHatAll, ySeqAll;
    vec lySamples;
    mat nuSamples, ltSamples;

    //Outsample variable
    mat dataXSeq, dataXRealTrans, dataXSeqTrans;
    int seqN;
    bool outSampleRequired;
    mat thetaSeq;
    vec ySeq;

    // prior parameters, the same for all component
    double ay, by, atheta, btheta, mu_theta;
    bool hasPrior;

    // related internal functions
    void oneRoundSample(bool recordStage);
    void adjustStepSize(int mcmcI);
    void initRecordSample(int nSample);
    void recordSample(int mcmcI);

    // Read the documentation for exp(S1)
    // Y_ is the model prediction value
    double compute_expS1(const vec &Y_)
    {
        double expS1_result;
        vec yDiff;
        yDiff = dataY - Y_;
        expS1_result = -0.5 * sum(yDiff % yDiff);
        return expS1_result;
    }

    // Read the documentation for exp(S2)
    double compute_expS2(const vec &theta_, const mat &U_, const vec &Lambda_)
    {
        double expS2_result;
        vec thetaDiff;
        thetaDiff = theta_ - mu_theta;
        thetaDiff = U_.t() * thetaDiff;
        thetaDiff /= sqrt(Lambda_);
        expS2_result = -0.5 * sum(thetaDiff % thetaDiff);
        return expS2_result;
    }

    //component 1, expS = exponential shoulder
    int Count1, Count2;
    double proposeC1;
    vec nu;
    vector<mat> U;
    vector<vec> Lambda;
    void sample_nu(int compI);
    double computeLikS_nu(double nu_, const mat &U_,
                          const vec Lambda_, const int compI);
    mat compute_kernel(double nu_);

    // component 2
    double proposeC2;
    vec yHat;
    mat thetaVec, thetaVecNew;
    void sample_thetaFun(int compI);
    vec compute_yHat(const mat &theta_);
    mat propose_thetaFun(int compI);

    // Component 3/4
    double lambda_y;
    vec lambda_theta;
    void sample_lambda_y();
    void sample_lambda_theta(int compI);

    // Commponent 5
    void compute_ySeq();
    void compute_thetaSeq(int compI);
};

#endif
