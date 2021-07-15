rm(list = ls())
source("./code/core/loadAll.R")
source("./code/simu1-parameter.R")


c1 = 1 # Simu 1, Cheap code
sfInit(parallel = TRUE, cpus = 16)
sfLibrary(calibration)
sfExportAll()
sfClusterSetupRNG(type="RNGstream", seed = 1231)
Simu1_1 = sfClusterApplyLB(1:100, oneRound_noCV)
sfStop()
save(Simu1_1,  file = "./data/simu1_1.RData")

