rm(list = ls())
#.rs.restartR()
library(calibration)
library(snowfall)
source("./code/core/multi_loadAll.R")
# compute MCMC evaluate loss
problemIndex = 1
thetaIndex = 6
linkIndex = 1
numComp = 2
linkLowerB = -20
linkUpperB = 20
lowerB = 1
upperB = 2
emuData = NULL
emulator_betaOpt = NULL
y_physics = y_physics_simu4
y_emulator = y_emulator_simu4
TrueSeqX = seq(lowerB, upperB, length.out = 200)
TrueSeqX = matrix(TrueSeqX, ncol = 1)
TrueSeqY = 1 + TrueSeqX^3
TrueSeqTheta = cbind(1/TrueSeqX, TrueSeqX)
nugget = 100

lenLambda = 50
lambdaSeq = rev(exp(seq(-7,3, length.out = lenLambda)))
GCVSeq = rep(0, lenLambda)
nSample = 50
nTry_global1 = 30
nTry_global2 = 30


mcmc_execute = function(calObj){
    calObj$setPrior(5, 5, 5, 5, 0)
    calObj$init_mcmc(0.5, 0.5, 1, 1, 1)
    res = calObj$execute_mcmc(3000, 3000)
    return(res)
}


### Simulation 3, Cheap Code
sfInit(parallel = TRUE, cpus = 16)
sfExportAll()
sfLibrary(calibration)
resultSF = sfClusterApplyLB(1:100,  oneRoundWrap_multi)
sfStop()
save(resultSF, file = "./data/simu3-1.RData")


