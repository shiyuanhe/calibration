#include "methodBase.hpp"


methodBase::methodBase(int c1, int c2, int c3, int numComp_){
    numComp = numComp_;
    sampleSize = 0;
    // select the objective function for the univariate case
    if(numComp == 1)
        selectSimulator(c1);
    
    thetaCoreVec.resize(numComp);
    thetaLinkVec.resize(numComp);
    size_t i;
    for (i = 0; i < numComp; i++)
    {
        selectThetaType(i, c2);
    }
    
    for (i = 0; i < numComp; i++)
    {
        selectLinkType(i, c3);
    }
    
}



void methodBase::selectSimulator(int c1){
    switch(c1){
    case 1:
        simulator = std::make_shared<simuobj_simu1>(); 
        break;
    case 2:
        simulator = std::make_shared<simuobj_simu2>(); 
        break;
    default:
        simulator = std::make_shared<simuobj_GP>(); break;
    }
}


void methodBase::selectThetaType(int i, int c2){
    switch (c2)
    {
    case 1:
        thetaCoreVec[i] = std::make_shared<thetaConst>();
        break;
    case 2:
        thetaCoreVec[i] = std::make_shared<thetaExponential>();
        break;
    case 3:
        thetaCoreVec[i] = std::make_shared<thetaQuadratic>();
        break;
    case 4:
        thetaCoreVec[i] = std::make_shared<thetaRKSH_linear_v1>();
        break;
    case 5:
        thetaCoreVec[i] = std::make_shared<thetaRKSH_cubic_v1>();
        break;
    }
}


void methodBase::selectLinkType(int i, int c3){
    switch (c3)
    {
    case 1:
        thetaLinkVec[i] = std::make_shared<linkIdentity>();
        break;
    case 2:
        thetaLinkVec[i] = std::make_shared<linkLogit>();
        break;
    case 3:
        thetaLinkVec[i] = std::make_shared<expLink>();
        break;
    }
}


void methodBase::set_linkBound(double lowerB_, double upperB_){
    size_t i;
    for(i = 0; i< numComp; i++)
        thetaLinkVec[i]->set_bound(lowerB_, upperB_);
}




// observed data, dataX has observed lower and upper bound
void methodBase::setData(vec dataY_, mat dataX_, vec lowerB_, vec upperB_){
    sampleSize = dataX_.n_rows;
    dataY = dataY_;
    dataXReal = dataX_;
    size_t i;
    
    for (i = 0; i < numComp; i++)
    {
        thetaCoreVec[i]->set_data(dataX_, lowerB_, upperB_);
    }
    
}
