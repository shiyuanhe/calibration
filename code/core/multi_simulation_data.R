## Multivariate true and emulator
y_physics_simu4 = function(nSample, lowerB, upperB){
    xData = matrix(runif(nSample, lowerB, upperB), ncol = 1)
    yPhysical = 1 + xData^3 + rnorm(nSample, sd = 0.2)
    return(list(x = xData, y = yPhysical))
}

y_emulator_simu4 = function(lowerB, upperB,linkLowerB, linkUpperB){
    thetaSeq = seq(linkLowerB, linkUpperB, length.out = 10)
    expGrid = expand.grid(thetaSeq, thetaSeq)
    expGrid = as.matrix(expGrid)
    
    xSeq = seq(lowerB, upperB, length.out = 10)
    expGrid = merge(xSeq, expGrid)
    expGrid = as.matrix(expGrid)
    
    ySimuGrid = expGrid[,2] * expGrid[,1] + expGrid[,3] * expGrid[,1]^2
    return(list(x = expGrid, y = ySimuGrid))
}



y_physics_simu5 = function(nSample, lowerB, upperB){
    xData = matrix(runif(nSample, lowerB, upperB), ncol = 1)
    yPhysical = xData^3 + rnorm(nSample, sd = 0.2)
    return(list(x = xData, y = yPhysical))
}


y_emulator_simu5 = function(lowerB, upperB,linkLowerB, linkUpperB){
    thetaSeq = seq(linkLowerB, linkUpperB, length.out = 10)
    expGrid = expand.grid(thetaSeq, thetaSeq)
    expGrid = as.matrix(expGrid)
    
    xSeq = seq(lowerB, upperB, length.out = 10)
    expGrid = merge(xSeq, expGrid)
    expGrid = as.matrix(expGrid)
    
    ySimuGrid = expGrid[,2]*expGrid[,1]^expGrid[,3]
    return(list(x = expGrid, y = ySimuGrid))
}

