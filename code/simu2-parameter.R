

thetaTrueFun = function(xx){
     0.5*(xx-2)^2 + 0.5
}
yTrueFun = function(xx){
     cos(2*xx) * sin(0.5*xx)
}



mcmc_execute = function(calObj){
    calObj$setPrior(8, 8, 8, 2, 0)
    #thetaInit = rep(-1, nSample)
    #calObj$init_mcmc(0.1, 0.1, 0, 1, 1, thetaInit)
    calObj$init_mcmc(0.1, 0.3, 1, 100, 2)
    testS = calObj$execute_mcmc(6e3, 6e3)
    return(testS)
}

## setting for simulation 2-1
c1 = 2
lowerB = 0.5 * pi
upperB = 1 * pi
nSample = 50

y_physics = y_physics_simu2
y_emulator = y_emulator_simu2
emulator_betaOpt = c(7.203282,-3.409236)
nugget = 1
lambdaSeq = rev(exp(seq(-15,-6,length.out = 30)))
lambdaSeqCubic = rev(exp(seq(-18,-12,length.out = 20)))


TrueSeqX = seq(lowerB, upperB, length.out = 200)
TrueSeqX = matrix(TrueSeqX, nrow = 200)
TrueSeqTheta = thetaTrueFun(TrueSeqX)
TrueSeqY = yTrueFun(TrueSeqX)
linkLowerB = min(TrueSeqTheta) * 0.7
linkUpperB =  max(TrueSeqTheta) * 1.1


## Overwrite Sampling Init

getInitGamma_cal1 = function(nSample){
    gammaInit = rnorm(1, mean = 0, sd = 2)
    return(gammaInit)
}


getInitGamma_cal2 = function(nSample){
    thetaValue = runif(2, linkLowerB, linkUpperB)
    x1 = 0; x2 = 1
    gamma1 = - log(abs(thetaValue[2]/thetaValue[1])) / 
        (x2 - x1)
    gamma0 = thetaValue[1] * exp(gamma1 * x1)
    return(c(gamma0, gamma1))
}

getInitGamma_cal3 = function(nSample){
    if(runif(1) < 0.1) return(rnorm(3, mean = c(2.5, -5.5, 5.5), sd = 0.5))
    thetaValue = sort(runif(3, linkLowerB, linkUpperB))
    x1 = 0; x3 = 1; x2 = (x1 + x3)/2
    A = rbind(c(1, x1, x1^2), c(1, x2, x2^2), c(1, x3, x3^2))
    gammaInit = solve(A, thetaValue)
    return(gammaInit)
}


getInitGamma_cal4 = function(nSample){
    gammaInit = rnorm(nSample + 1, mean = 0, sd = 4)
    #gammaInit[1] = -1
    return(gammaInit)
}
# 
getInitGamma_cal6 = function(nSample){
    gammaInit = rnorm(nSample + 2, mean = 0, sd = 5)
     gammaInit[1] = -5 + gammaInit[1] 
     gammaInit[2] = -3 + gammaInit[2]
    return(gammaInit)
}


