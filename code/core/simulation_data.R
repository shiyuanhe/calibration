

y_physics_simu1 = function(nSample, lowerB, upperB){
    xData = runif(nSample, lowerB, upperB)
    yData = exp(xData / 10) * cos(xData)  + rnorm(nSample, sd = 0.1)
    xData = matrix(xData, nrow = nSample)
    yData = matrix(yData, nrow = nSample)
    return(list(x = xData, y = yData))
}

y_emulator_simu1 = function(lowerB, upperB, linkLowerB, linkUpperB){
    xSeq = seq(lowerB, upperB, length.out = 14)
    thetaSeq = seq(linkLowerB, linkUpperB, length.out = 15)
    expGrid = expand.grid(xSeq, thetaSeq)
    expGrid = as.matrix(expGrid)
    ySimuGrid = exp(expGrid[,1]/10) * cos(expGrid[,1]) * exp(expGrid[,1]/5)/expGrid[,2] * 0.5
    return(list(x = expGrid, y = ySimuGrid))
}




y_physics_simu2 = function(nSample, lowerB, upperB){
    xData = runif(nSample, lowerB, upperB)
    yData =  cos(2*xData) * sin(0.5 * xData)  + rnorm(nSample, sd = 0.1)
    xData = matrix(xData, nrow = nSample)
    yData = matrix(yData, nrow = nSample)
    return(list(x = xData, y = yData))
}

y_emulator_simu2 = function(lowerB, upperB,linkLowerB, linkUpperB){
    xSeq = seq(lowerB, upperB, length.out = 14)
    thetaSeq = seq(linkLowerB , linkUpperB , length.out = 15)
    expGrid = expand.grid(xSeq, thetaSeq)
    expGrid = as.matrix(expGrid)
    
    yData = cos(2*expGrid[,1]) * sin(0.5 * expGrid[,1])
    den = 0.5 *(expGrid[,1] - 2)^2 + 0.5
    ratio = expGrid[,2] / den
    #ratio[ratio>2] = 2
    #ratio[ratio<(-2)] = -2
    yData = yData * sin(0.5*pi) * cos(0*pi) * exp(3*ratio-3)
    
    return(list(x = expGrid, y = yData))
}

