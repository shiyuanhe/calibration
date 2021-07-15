computeErr = function(yHat_OOB){
    errConst = mean(abs(yHat_OOB - real_yp))
    errConstSD = sd(abs(yHat_OOB - real_yp))
    c(errConst, errConstSD/sqrt(nSample))
}



getInitGamma_cal6 = function(nSample){
    gammaInit = rnorm(nSample + 2, mean = 0, sd = 10)
    return(gammaInit)
}



RKHS_OOB = function(c2, lambdaSeq, getInitGamma, c3 = 2){
    
    yHat_OOB = rep(0, nSample)
    for(rmI in 1:nSample){
        # cat(rmI, "\r")
        # flush.console()
        r1Obj = create_CalibrationObj(problemIndex = 4, thetaIndex = c2, linkIndex = c3, 
                                      numComp = 1, 
                                      emuData = emuData, betaOpt = emu_betaOpt, 
                                      linkLowerB = linkLowerB, linkUpperB = linkUpperB)
        real_x_cv = matrix(real_x[-rmI,], nrow = nSample - 1)
        r1Obj$setData(real_yp[-rmI], real_x_cv,  lowerB, upperB)
        result = selectLambdaGCV(r1Obj, lambdaSeq, 
                                 nSample - 1, 
                                 getInitGamma, TrueSeqX, 
                                 linkLowerB, linkUpperB,
                                 plotcv = FALSE)
        cv_x = matrix(real_x[rmI,], nrow = 1)
        yHat_OOB[rmI] = r1Obj$predict_y(result$gammaOpt, cv_x, FALSE)[1,1]
    }
    return(yHat_OOB)
}


# leave two out cross-validation
RKHS_OOB_randOUT = function(c2, nRep, lambdaSeq, getInitGamma, c3 = 2, oobSize){
    #oobSize = 2#floor(nSample * 0.3)
    trainSize = nSample - oobSize
    yHat_OOB = matrix(0, nRep,1)
    for(repI in 1:nRep){
        # cat(rmI, "\r")
        # flush.console()
        rmSet = sample(1:nSample, oobSize,replace = FALSE)
        r1Obj = create_CalibrationObj(problemIndex = 4, thetaIndex = c2, linkIndex = c3,
                                      numComp = 1,
                                      emuData = emuData, betaOpt = emu_betaOpt,
                                      linkLowerB = linkLowerB, linkUpperB = linkUpperB)
        real_x_train = real_x[-rmSet,,drop = FALSE]
        r1Obj$setData(real_yp[-rmSet], real_x_train,  lowerB, upperB)
        result = selectLambdaGCV(r1Obj, lambdaSeq,
                                 trainSize,
                                 getInitGamma, TrueSeqX,
                                 linkLowerB, linkUpperB,
                                 plotcv = FALSE)
        cv_x = matrix(real_x[rmSet,], nrow = oobSize)
        resids = r1Obj$predict_y(result$gammaOpt, cv_x, FALSE)[,1] - real_yp[rmSet]
        yHat_OOB[repI,1] = sum(abs(resids))
     
    }
    return(yHat_OOB)
}


RKHS_RD = function(c2, lambdaSeq, getInitGamma, predictX, c3 = 2){

    
    r1Obj = create_CalibrationObj(problemIndex = 4, thetaIndex = c2, linkIndex = c3, 
                                  numComp = 1, 
                                  emuData = emuData, betaOpt = emu_betaOpt, 
                                  linkLowerB = linkLowerB, linkUpperB = linkUpperB)
    r1Obj$setData(real_yp, real_x, lowerB, upperB)
    result = selectLambdaGCV(r1Obj, lambdaSeq, 
                             nSample, 
                             getInitGamma, TrueSeqX, 
                             linkLowerB, linkUpperB,
                             plotcv = TRUE)
    fitX = r1Obj$predict_y(result$gammaOpt, real_x, TRUE)
    fitNew = r1Obj$predict_y(result$gammaOpt, predictX, TRUE)
    
    return(list(fitNew, fitX))
}



plotData1 = function(thetaSeq, spath = NULL){
    if(!is.null(spath)){
        pdf(spath, width = 4, height = 3)
    }
    par(mar = c(4,4,1,1))
    plot(predX, thetaSeq[,3], type = "l", 
         xlab = "PVA Amount", ylab = "Effectiveness",
         ylim = c(0.65, 1.4))
    lines(predX, thetaSeq[,3] - thetaSeq[,4], type = "l", lty = 2)
    lines(predX, thetaSeq[,3] + thetaSeq[,4], type = "l", lty = 2)
    if(!is.null(spath)){
        dev.off()
    }
}


