rm(list = ls())
#.rs.restartR()
library(calibration)
library(snowfall)
source("./code/core/multi_loadAll.R")
# compute MCMC evaluate loss
problemIndex = 2
thetaIndex = 5
linkIndex = 1
numComp = 2
linkLowerB = 0
linkUpperB = 4
lowerB = 1
upperB = 2
emuData = NULL
emulator_betaOpt = NULL
y_physics = y_physics_simu5
y_emulator = y_emulator_simu5
TrueSeqX = seq(lowerB, upperB, length.out = 200)
TrueSeqX = matrix(TrueSeqX, ncol = 1)
TrueSeqY = TrueSeqX^3
TrueSeqTheta = cbind(1, TrueSeqX)
nugget = 100

lenLambda = 50
lambdaSeq = rev(exp(seq(-5,5, length.out = lenLambda)))
GCVSeq = rep(0, lenLambda)
nSample = 50
nTry_global1 = 30
nTry_global2 = 30


mcmc_execute = function(calObj){
    # a_y, b_y, a_theta, b_theta, mu_teta
    calObj$setPrior(5, 5, 5, 5, 0)
    # stepsize
    calObj$init_mcmc(0.5, 0.5, 1, 1, 1)
    res = calObj$execute_mcmc(3000,3000)
    return(res)
}

sfInit(parallel = TRUE, cpus = 16)
sfExportAll()
sfLibrary(calibration)
resultSF2 = sfClusterApplyLB(1:100,  oneRoundWrap_multi)
sfStop()

save(resultSF2, file = "./data/simu4-1.RData")


