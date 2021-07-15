## Get random initial values for theta

getInitGamma_cal1 = function(nSample){
    gammaInit = rnorm(1, mean = 0, sd = 2)
    return(gammaInit)
}

getInitGamma_cal2 = function(nSample){
    if(runif(1) < 0.1) return(c(1/2, -1/5))
    gammaInit = rnorm(2, mean = 0, sd = c(1, 0.5))
    return(gammaInit)
}

getInitGamma_cal3 = function(nSample){
    if(runif(1) < 0.1) return(c(1, -4, 5))
    gammaInit = rnorm(3, mean = 0, sd = 3)
    return(gammaInit)
}

getInitGamma_cal4 = function(nSample){
    gammaInit = rnorm(nSample + 1, mean = 0, sd = 10)
    return(gammaInit)
}

getInitGamma_cal6 = function(nSample){
    gammaInit = rnorm(nSample + 2, mean = 0, sd = 10)
    gammaInit[1:2] = -1
    return(gammaInit)
}

getInitGamma_cal8 = function(nSample){
    gammaInit = rnorm(nSample, mean = 0, sd = 10)
    gammaInit[1:2] = -1
    return(gammaInit)
}