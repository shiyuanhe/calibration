## One round of simulation for the Bayesian method
## Many things should reside in the global memeory before computing
##     TrueSeqX, TrueSeqTheta, TrueSeqY
##     mcmc_execute()

alphaList_mcmc = c(0.05, 0.025, 0.005)
oneRound_bayesian = function(iterI){
    emuData = y_emulator(lowerB, upperB, linkLowerB, linkUpperB)
    simuData = y_physics(nSample,lowerB, upperB)
    
    calObjMCMC = mcmc_CalibrationObj(problemIndex, thetaIndex, linkIndex,
                                     numComp, linkLowerB, linkUpperB,
                                     emuData, emulator_betaOpt)
    calObjMCMC$setData(simuData$y, simuData$x, lowerB, upperB)
    calObjMCMC$setPredition(TrueSeqX)
    mcmcResult = mcmc_execute(calObjMCMC)
    #ThetaHat = apply(mcmcResult$thetaSeqAll, c(1,2), mean)
    loss1 = bayesianLoss(mcmcResult)
    # res = c(99, loss1)
    loss2 = bayesianPredictionLoss(mcmcResult)
    #loss2 = c(99, loss2)
    
    result = rbind(loss2,loss1)
    ntheta = nrow(loss1)
    rownames(result) = c("Y",paste0("theta",1:ntheta))
    tmp_colnames = outer(c("Width","Coverage") ,1-2*alphaList_mcmc, FUN = "paste0")
    dim(tmp_colnames) = NULL
    colnames(result) = c("loss",tmp_colnames)
    return(result)
}


## The Bayesian Prediction Error
bayesianPredictionLoss = function(mcmcResult){
    # yHat = rowMeans(mcmcResult$ySeqAll)
    # error = sum((yHat - TrueSeqY)^2) * diff(TrueSeqX)[1]
    # return(error)

    applyDim = 1#c(1,2)
    yHat = apply(mcmcResult$ySeqAll, applyDim, mean)
    deltaX = diff(TrueSeqX)[1]
    diff = TrueSeqY - yHat
    L2Loss = sqrt(sum(diff^2) * deltaX)

    res = L2Loss

    ## Compute the upper and lower bound, then compute error
    for(alpha in alphaList_mcmc){
        UpperB = apply(mcmcResult$ySeqAll, applyDim, quantile, probs = 1 - alpha)
        LowerB = apply(mcmcResult$ySeqAll, applyDim, quantile, probs = alpha)
        coverB = LowerB < TrueSeqY & TrueSeqY < UpperB
        cm1 = sum(coverB) / length(coverB)
        #cm2 = all(coverB)
        cWidth = sum( (UpperB - LowerB)) * deltaX
        res = c(res, cWidth, cm1)
    }
    return(res)
}

# For the computed result of the Bayesian method, 
# compute its L2 error and coverage error. 
bayesianLoss = function(mcmcResult){
    
    applyDim = c(1,2)
    ThetaHat = apply(mcmcResult$thetaSeqAll, c(1,2), mean)
    deltaX = diff(TrueSeqX)[1]
    diff = TrueSeqTheta - ThetaHat
    L2Loss = sqrt(colSums(diff^2) * deltaX)
    
    res = L2Loss
    
    ## Compute the upper and lower bound, then compute error
    tmpRes = L2Loss
    for(alpha in alphaList_mcmc){
        UpperB = apply(mcmcResult$thetaSeqAll, applyDim, quantile, probs = 1 - alpha)
        LowerB = apply(mcmcResult$thetaSeqAll, applyDim, quantile, probs = alpha)
        tmp = coverageMeasure(LowerB, UpperB, TrueSeqTheta)
        tmpRes = cbind(tmpRes,tmp)
    }
    return(tmpRes)
}

## Compute the width of the confidence reigon, coverage proportion, 
## whether or not it is fully coveraged. 
coverageMeasure = function(LowerB, UpperB, TrueSeqTheta){
    deltaX = diff(TrueSeqX)[1]
    diffB = UpperB - LowerB
    CIWidth = colSums(diffB) * deltaX
    
    coverB = LowerB < TrueSeqTheta & TrueSeqTheta < UpperB
    cm1 = colSums(coverB) / nrow(coverB)
    cm2 = apply(coverB, 2, all)
    return(cbind(CIWidth, cm1)) #,cm2
}



## Check if the whole curve is covered
# bayesianPredictionLossWhole = function(mcmcResult){
#     # yHat = rowMeans(mcmcResult$ySeqAll)
#     # error = sum((yHat - TrueSeqY)^2) * diff(TrueSeqX)[1]
#     # return(error)
#     
#     applyDim = 1#c(1,2)
#     yHat = apply(mcmcResult$ySeqAll, applyDim, mean)
#     deltaX = diff(TrueSeqX)[1]
#     diff = TrueSeqY - yHat
#     L2Loss = sqrt(sum(diff^2) * deltaX)
#     
#     res = L2Loss
#     
#     ## Compute the upper and lower bound, then compute error
#     for(alpha in alphaList_mcmc){
#         UpperB = apply(mcmcResult$ySeqAll, applyDim, quantile, probs = 1 - alpha)
#         LowerB = apply(mcmcResult$ySeqAll, applyDim, quantile, probs = alpha)
#         coverB = LowerB < TrueSeqY & TrueSeqY < UpperB
# 
#         
#         cm1 = sum(coverB) / length(coverB)
#         cm2 = all(coverB)
#         cWidth = sum( (UpperB - LowerB)) * deltaX
#         res = c(res, cWidth, cm2)
#     }
#     return(res)
# }

# Check if the whole curve is covered
bayesianLossWhole = function(mcmcResult){
    
    applyDim = c(1,2)
    ThetaHat = apply(mcmcResult$thetaSeqAll, c(1,2), mean)
    deltaX = diff(TrueSeqX)[1]
    diff = TrueSeqTheta - ThetaHat
    L2Loss = sqrt(colSums(diff^2) * deltaX)
    
    res = L2Loss
    
    ## Compute the upper and lower bound, then compute error
    tmpRes = L2Loss
    zNum = dim(mcmcResult$thetaSeqAll)[3]
    zVec = rep(0,zNum)
    for(i in 1:zNum){
        zVec[i] = max(abs(mcmcResult$thetaSeqAll[,,i] - ThetaHat))
    }
    
    for(alpha in alphaList_mcmc){
        # UpperB = apply(mcmcResult$thetaSeqAll, applyDim, quantile, probs = 1 - alpha)
        # LowerB = apply(mcmcResult$thetaSeqAll, applyDim, quantile, probs = alpha)
        shift = quantile(zVec,1-2*alpha)
        UpperB = ThetaHat + shift
        LowerB = ThetaHat - shift

        
        deltaX = diff(TrueSeqX)[1]
        diffB = UpperB - LowerB
        CIWidth = colSums(diffB) * deltaX
        
        coverB = LowerB < TrueSeqTheta & TrueSeqTheta < UpperB
        cm2 = apply(coverB, 2, all)
        
        tmp = cbind(CIWidth, cm2)
        
        
        tmpRes = cbind(tmpRes,tmp)
    }
    return(tmpRes)
}

