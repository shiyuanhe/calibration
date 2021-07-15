## Get random initial values for theta


thetaIndexList = c(1,2,3,6)
getInit1 = function(nSample){
    gamma = matrix(rnorm(2, sd = 5), ncol = 2)
    return(gamma)
}
getInit2 = function(nSample){
    gamma = matrix(rnorm(4, sd = 5), ncol = 2)
    return(gamma)
}
getInit3 = function(nSample){
    gamma = matrix(rnorm(6, sd = 5), ncol = 2)
    return(gamma)
}
getInit6 = function(nSample){
    gamma = matrix(rnorm(2 * (nSample + 2), sd = 5), ncol = 2)
    return(gamma)
}

getInitList = list(getInit1, getInit2, getInit3, getInit6)



