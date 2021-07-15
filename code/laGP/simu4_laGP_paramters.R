# The laGP parameters for the third simulation in the paper



# cheap code
cheapModel = function(X, THETA){
    X <- as.matrix(X)
    THETA <- as.matrix(THETA)
    
    #####ySimuGrid = expGrid[,2] * expGrid[,1] + expGrid[,3] * expGrid[,1]^2
    #Y<- THETA[,1]*X[,1] + THETA[,2]*X[,1]^2
    
    ### ySimuGrid = expGrid[,2]*expGrid[,1]^expGrid[,3]
    Y <- THETA[,1]*X[,1]^THETA[,2]
    return(Y)
}


beta.prior = NULL
nThread = 6
## set up priors
bias.est <- TRUE # estimate the bias 
methods <- rep("alc", 2)

xDim = 1 # dimension of the predictor x
snomadr_n = ncalib = 2 # dimension of the calibration parameter
snomadr_bbin = c(0,0) # variable type (0: continous) for the input 
snomadr_bbout = 0 # variable type (0: continous) for the output 
snomadr_lb = c(0,0) # lower bound for the calibration parameter
snomadr_ub = c(4,4) # uppder bound for the calibration paramter

imesh = 0.1
opts <- list("MAX_BB_EVAL"=1000, "INITIAL_MESH_SIZE"=imesh, "MIN_POLL_SIZE"="r0.001")


hyperD = 0.5
hyperG = NULL