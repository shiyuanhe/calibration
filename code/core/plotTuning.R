plotBoot = function(TrueSeqX, TrueSeqTheta, sList, thetaHat){
    plot(TrueSeqX, TrueSeqTheta, type = "n", col = "blue", 
         ylim = c(linkLowerB, linkUpperB))
    upperB = sList$U90
    lowerB = sList$L90
    polygon(c(TrueSeqX, rev(TrueSeqX)), c(upperB, rev(lowerB)),
            col = adjustcolor("red", alpha.f=0.3), lty = 0)
    lines(TrueSeqX, TrueSeqTheta, col = "blue", lwd = 2)
    lines(TrueSeqX, thetaHat,  col = "red", lwd = 2)
    #lines(TrueSeqX, sList$thetaHat, lty = 2, col = "red")
    
}


yHat2CBList =  function(yHat){
    upperB90 = yHat[,1] + yHat[,2] * qnorm(0.95)
    lowerB90 = yHat[,1] - yHat[,2] * qnorm(0.95)
    upperB95 = yHat[,1] + yHat[,2] * qnorm(0.975)
    lowerB95 = yHat[,1] - yHat[,2] * qnorm(0.975)
    upperB99 = yHat[,1] + yHat[,2] * qnorm(0.995)
    lowerB99 = yHat[,1] - yHat[,2] * qnorm(0.995)
    return(CBList = list(U90 = upperB90, U95 = upperB95, U99 = upperB99, 
                         L90 = lowerB90, L95 = lowerB95, L99 = lowerB99))
}

plotModelFitting = function(yHat, ySigma, thetaHat, CBList, simuData,
                            thetaHat2 = NULL){
    par0 = par(mfrow=  c(2,1), mar = c(2,2,1,1))
    plot(TrueSeqX, TrueSeqTheta, type = "n", col = "blue", 
         ylim = c(linkLowerB, linkUpperB))
    
    polygon(c(TrueSeqX, rev(TrueSeqX)), c(CBList$U90, rev(CBList$L90)),
            col = adjustcolor("green", alpha.f=0.3), lty = 0)
    polygon(c(TrueSeqX, rev(TrueSeqX)), c(CBList$U95, rev(CBList$U90)),
            col = adjustcolor("blue", alpha.f=0.6), lty = 0)
    polygon(c(TrueSeqX, rev(TrueSeqX)), c(CBList$U95, rev(CBList$U99)),
            col = adjustcolor("red", alpha.f=0.8), lty = 0)
    polygon(c(TrueSeqX, rev(TrueSeqX)), c(CBList$L95, rev(CBList$L90)),
            col = adjustcolor("blue", alpha.f=0.6), lty = 0)
    polygon(c(TrueSeqX, rev(TrueSeqX)), c(CBList$L99, rev(CBList$L95)),
            col = adjustcolor("red", alpha.f=0.8), lty = 0)
    
    lines(TrueSeqX, thetaHat, pch = 20, col = "red", lwd = 1)
    lines(TrueSeqX, TrueSeqTheta, col = "black", lwd = 3)
    if(!is.null(thetaHat2)){
        lines(TrueSeqX, thetaHat2, col = "red", lwd = 1, lty = 2)
    }
    ## yHat and its uncertainty
    upperB = yHat + ySigma * 1.645
    lowerB = yHat - ySigma * 1.645
    plot(simuData$x, simuData$y, pch = 20, col = "blue")
    polygon(c(TrueSeqX, rev(TrueSeqX)), c(upperB, rev(lowerB)),
            col = adjustcolor("red", alpha.f=0.3), lty = 0)
    lines(TrueSeqX, yHat, pch = 20, col = "red", lwd = 2)
    lines(TrueSeqX, TrueSeqY, type = "l", col = "blue", lwd = 2)
    par(par0)
}