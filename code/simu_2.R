rm(list = ls())
source("./code/core/loadAll.R")
source("./code/simu2-parameter.R")


c1 = 2 # Simulation 2 as objective model
sfInit(parallel = TRUE, cpus = 16)
sfClusterSetupRNG(type="RNGstream", seed = 1231)
sfLibrary(calibration)
sfExportAll()
Simu2_1 = sfClusterApplyLB(1:100, oneRound_noCV)
sfStop()
save(Simu2_1,  file = "./data/simu2_1.RData")
