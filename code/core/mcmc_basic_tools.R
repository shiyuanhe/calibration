## create mcmc object
mcmc_CalibrationObj = function(problemIndex, thetaIndex, linkIndex, numComp,
                               linkLowerB, linkUpperB,
                               emuData = NULL, betaOpt = NULL){
    calObj = new(mcmcCal, problemIndex, thetaIndex, linkIndex, numComp)
    calObj$set_linkBound(linkLowerB, linkUpperB)
    if(!is.null(emuData) & !is.null(betaOpt)){
        calObj$GP_setData(emuData$y, emuData$x)
        calObj$GP_setBeta(betaOpt)
    }
    return(calObj)
}


mcmc_predErr = function(calObj, lowerB, upperB,
                           trainingData, testingData){
    calObj$setData(trainingData$y, trainingData$x,  lowerB, upperB)
    resultMCMC = mcmc_execute(calObj)
    ySeq = rowMeans(resultMCMC$ySeqAll)
    yFun = approxfun(TrueSeqX, ySeq)
    yHat = yFun(testingData$x)
    predErr = yHat - testingData$y
    return(predErr)
}

mcmc_thetaErr = function(calObj, lowerB, upperB,
                        trainingData){
    calObj$setData(trainingData$y, trainingData$x,  lowerB, upperB)
    resultMCMC = mcmc_execute(calObj)
    thetaHatSeq = matrixMeanFunction(resultMCMC$thetaSeqAll)
    Err1 = TrueSeqTheta - thetaHatSeq
    Err1 = sum(Err1^2) * diff(TrueSeqX)[1]
    Err1 = sqrt(Err1)
    
    yHatSeq = matrixMeanFunction(resultMCMC$ySeqAll) ## ??? Error
    Err2 = TrueSeqY - yHatSeq
    Err2 = sum(Err2^2) * diff(TrueSeqX)[1]
    Err2 = sqrt(Err2)
    return(c(Err1, Err2))
}

mcmc_thetaErr_empirical = function(calObj, lowerB, upperB,
                         trainingData){
    calObj$setData(trainingData$y, trainingData$x,  lowerB, upperB)
    resultMCMC = mcmc_execute(calObj)
    thetaHatSeq = matrixMeanFunction(resultMCMC$thetaS)
    Err1 = thetaTrueFun(trainingData$x) - thetaHatSeq
    Err1 = sqrt(mean(Err1^2))

    yHatSeq = rowMeans(resultMCMC$yHat)
    Err2 = yTrueFun(trainingData$x) - yHatSeq
    Err2 = sqrt(mean(Err2^2))
    
    return(c(Err1, Err2))
}

