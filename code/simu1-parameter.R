
## setting for simulation 1-1
thetaTrueFun = function(xx){
    exp(xx/5) / 2
}
yTrueFun = function(xx){
    exp(xx / 10) * cos(xx)
}

mcmc_execute = function(calObj){
    calObj$setPrior(8, 8, 2, 1, 0)
    calObj$init_mcmc(0.1, 0.3, 1, 100, 2)
    testS = calObj$execute_mcmc(6e3, 6e3)
    return(testS)
}



##
c1 = 1
lowerB = 1 * pi
upperB = 3 * pi
nSample = 50
emulator_betaOpt = c(1.1126047, -0.9420013)
nugget = 1
lambdaSeq = rev(exp(seq(-16,-6,length.out = 13)))
lambdaSeqCubic = (exp(seq(-16, -12,length.out = 30)))

##
y_physics = y_physics_simu1
y_emulator = y_emulator_simu1
TrueSeqX = seq(lowerB, upperB, length.out = 200)
TrueSeqX = matrix(TrueSeqX, nrow = 200)
TrueSeqTheta = thetaTrueFun(TrueSeqX)
TrueSeqY = yTrueFun(TrueSeqX)
linkLowerB = min(TrueSeqTheta) * 0.7
linkUpperB =  max(TrueSeqTheta) * 1.1


##

## Get random initial values for theta

getInitGamma_cal1 = function(nSample){
    gammaInit = rnorm(1, mean = 0, sd = 2)
    return(gammaInit)
}


## Initialize around the true optimal once every 10 times
getInitGamma_cal2 = function(nSample){
    if(runif(1) < 0.1) return(rnorm(2, mean = c(1, -1), sd = 0.1))
    thetaValue = runif(2, linkLowerB, linkUpperB)
    x1 = 0; x2 = 1
    gamma1 = - log(abs(thetaValue[2]/thetaValue[1])) / 
        (x2 - x1)
    gamma0 = thetaValue[1] * exp(gamma1 * x1)
    return(c(gamma0, gamma1))
}

getInitGamma_cal3 = function(nSample){
    if(runif(1) < 0.1) return(rnorm(3, mean = 1, sd = 0.1))
    thetaValue = sort(runif(3, linkLowerB, linkUpperB))
    x1 = 0; x3 = 1; x2 = (x1 + x3)/2
    A = rbind(c(1, x1, x1^2), c(1, x2, x2^2), c(1, x3, x3^2))
    gammaInit = solve(A, thetaValue)
    return(gammaInit)
}


getInitGamma_cal6 = function(nSample){
    gammaInit = rnorm(nSample + 2, mean = 0, sd = 10)
    return(gammaInit)
}


