## To compare different methods
## Evaluate loss

## compute loss / posterior CI and bootstrap CI for
## calObjCC with tuning parameter lambda
## simuData$x in [lowerB, upperB]
## 200 TrueSeqX in [lowerB, upperB], with true theta values TrueSeqTheta
oneRoundCore = function(calObjCC, getInit, 
                        simuData, lowerB, upperB, 
                        TrueSeqX, TrueSeqTheta, nTry,
                        linkLowerB, linkUpperB){
    ## Initial estimate and confidence band from Posterior
    gammaOpt = optimize_CalModel(getInit, calObjCC, nSample, nTry,
                                 TrueSeqX, linkLowerB, linkUpperB)
    yHat = calObjCC$predict_y(gammaOpt, TrueSeqX, TRUE)
    ThetaHat = yHat[,3]
    ThetaSigma = yHat[,4]
    # CBList = yHat2CBList(yHat)
    # ComputeThetaLoss(TrueSeqX, TrueSeqTheta, yHat[,1], yHat[,2])
    loss = ComputeThetaLoss(TrueSeqX, TrueSeqTheta, ThetaHat, ThetaSigma)
    
    return(loss)
}



oneRound_noCV = function(iterI){
    c2 = c(1,2,3,6) # 6 cubic
    c3 = c(2,1,1,2) # 1 identity 2 logit bound
    getInitList = list(getInitGamma_cal1, getInitGamma_cal2, 
                       getInitGamma_cal3, getInitGamma_cal6)

    
    
    emuData = y_emulator(lowerB, upperB, linkLowerB, linkUpperB)
    simuData = y_physics(nSample,lowerB, upperB)
    
    predErrMat = numeric(0)
    for(j in 1:length(c2)){
        calObjCC = create_CalibrationObj(problemIndex = c1, 
                                         thetaIndex = c2[j], linkIndex = c3[j], 
                                         numComp = 1,
                                         emuData = emuData, 
                                         betaOpt = emulator_betaOpt, 
                                         linkLowerB, linkUpperB)
        calObjCC$setData(simuData$y, simuData$x,  lowerB, upperB)
        calObjCC$GP_setNugget(nugget)
        if(c2[j]==4)
            lambdaSeq = lambdaSeqLinear
        if(c2[j]==5)
            lambdaSeq = lambdaSeqCubic
        if(c2[j]==6)
            lambdaSeq = lambdaSeqCubic#lambdaSeqExp
        
        if(c2[j] > 3){
            result = selectLambdaGCV(calObjCC, lambdaSeq, 
                                     nSample, 
                                     getInitList[[j]], TrueSeqX, 
                                     linkLowerB, linkUpperB,
                                     plotcv = FALSE)
            calObjCC$setLambda(result$lambdaOpt)
            yHat = calObjCC$predict_y(result$gammaOpt, TrueSeqX, TRUE)


            res = ComputeThetaLoss(TrueSeqX, TrueSeqTheta, yHat[,3], yHat[,4]) #*sqrt(2)
            

        }else{
            nTry = c(50, 100, 100, 200, 200, 200)[j]
            res = oneRoundCore(calObjCC, getInitList[[j]], simuData,
                               lowerB, upperB,  TrueSeqX, 
                               TrueSeqTheta, nTry,
                               linkLowerB, linkUpperB)
        }
        predErrMat = rbind(predErrMat, res)

    }
    
    
 
    #compute MCMC evaluate loss
    calObjMCMC = mcmc_CalibrationObj(problemIndex = c1, thetaIndex =  1, linkIndex = 2,
                                     numComp = 1,emuData = emuData,
                                     betaOpt = emulator_betaOpt,
                                     linkLowerB = linkLowerB, linkUpperB = linkUpperB)
    calObjMCMC$setData(simuData$y, simuData$x, lowerB, upperB)
    calObjMCMC$setPredition(TrueSeqX)
    mcmcResult = mcmc_execute(calObjMCMC)
    loss1 = bayesianLoss( mcmcResult) #TrueSeqX, TrueSeqTheta,
     
   
    predErrMat = rbind(predErrMat, loss1[1,])
    

    
    return(predErrMat)
}
