#include "gaussianProcess.hpp"
#include "methodBayes.hpp"
#include "methodRKHS.hpp"
#include "theta_RKHS_cubic_v1.hpp"
#include "method_RKHS_multivariate.hpp"

RCPP_MODULE(computerCal)
{
    class_<methodBase>("methodBase")
        .constructor<int, int, int, int>()
        .method("setData", &methodBase::setData)
        .method("GP_setData", &methodBase::GP_setData)
        .method("GP_setBeta", &methodBase::GP_setBeta)
        .method("GP_setNugget", &methodBase::GP_setNugget)
        .method("GP_setScaling", &methodBase::GP_setScaling)
        .method("set_linkBound", &methodBase::set_linkBound)
        .method("getNumComp", &methodBase::getNumComp)
        //.method("setAuxData", &methodBase::setAuxData)
    ;

    class_<calModel>("calModel")
        .derives<methodBase>("methodBase")
        .constructor<int, int, int, int>()
        .method("setLambda", &calModel::setLambda)
        .method("objFun", &calModel::objFun)
        .method("gradFun", &calModel::gradFun)
        .method("predict_y", &calModel::predict_y)
        .method("fittedCovMat", &calModel::fittedCovMat)
        .method("compute_GCV", &calModel::compute_GCV)
        .method("compute_Sigma", &calModel::compute_residualSigma);

    class_<calModelMulti>("calModelMulti")
        .derives<methodBase>("methodBase")
        .constructor<int, int, int, int>()
        .method("setLambda", &calModelMulti::setLambda)
        .method("objFun", &calModelMulti::objFun)
        .method("gradFun", &calModelMulti::gradFun)
        .method("predict_y", &calModelMulti::predict_y)
        .method("sigmaSq_y", &calModelMulti::sigmaSq_y)
        .method("sigmaSq_y_simu", &calModelMulti::sigmaSq_y_simu)
        .method("compute_sigmaSq", &calModelMulti::compute_sigmaSq)
        .method("compute_GCV", &calModelMulti::compute_GCV);
    

    class_<mcmcCal>("mcmcCal")
        .derives<methodBase>("methodBase")
        .constructor<int, int, int, int>()
        .method("setPrior", &mcmcCal::setPrior)
        .method("execute_mcmc", &mcmcCal::execute_mcmc)
        .method("init_mcmc", &mcmcCal::init_mcmc)
        .method("setPredition", &mcmcCal::setPredition);

    class_<gaussianP>("gaussianP")
        .constructor()
        .method("set_beta", &gaussianP::set_beta)
        .method("cv_obj", &gaussianP::cv_obj)
        .method("cv_grad", &gaussianP::cv_grad)
        .method("predictY", &gaussianP::predict_y)
        .method("predictCov", &gaussianP::predict_yCov)
        .method("set_data", &gaussianP::set_data)
        .method("set_scaling", &gaussianP::set_scaling)
        .method("set_nugget", &gaussianP::set_nugget)
    ;

}
