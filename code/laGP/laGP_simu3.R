# Apply la GP to the third simulation
rm(list = ls())
#.rs.restartR()
library(calibration)
library(snowfall)
source("./code/core/multi_loadAll.R")
source("./code/laGP/laGP_core.R")
source("./code/laGP/simu3_laGP_paramters.R")

# compute MCMC evaluate loss
problemIndex = 1
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
nSample = 50




## cheap code where cheapModel is not NULL.
set.seed(100)
nRep = 100
resMatCheap = matrix(0,nRep,7)
for(i in 1:nRep){
    simuData = y_physics(nSample,lowerB, upperB)
    emuData = y_emulator(lowerB, upperB,linkLowerB, linkUpperB)
    res = laGP_predict(simuData$x,simuData$y, emuData,TrueSeqX)
    error = ComputeThetaLoss(TrueSeqX, TrueSeqY, res[,1], res[,2])
    resMatCheap[i,] = error
}

laGP_print(resMatCheap)


save(resMatCheap,file = "./data/simu3_laGP.RData")
