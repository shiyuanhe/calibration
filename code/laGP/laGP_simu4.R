# Apply la GP to the forth simulation in the manuscript
rm(list = ls())
#.rs.restartR()
library(calibration)
library(snowfall)
source("./code/core/multi_loadAll.R")
source("./code/laGP/laGP_core.R")
source("./code/laGP/simu4_laGP_paramters.R")

# compute MCMC evaluate loss
problemIndex = 2
thetaIndex = 6
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


save(resMatCheap,file = "./data/simu4_laGP.RData")

