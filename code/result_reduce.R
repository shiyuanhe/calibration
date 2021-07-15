# ----------- simu1 ------------
get_mean_value = function(SimuRes,  colT2,colT3){
    mValue = matrix(0, rowN, colN)
    sdValue = matrix(0, rowN, colN)
    mName = c("Const", "Param-Exp", "Param-Quad",
              "RKHS-Cubic", "Baysian")
    for(i in 1:rowN){
        # Computation
        for(j in colT2){
            tmpS = sapply(SimuRes, function(resM) resM[i,j])
            mValue[i,j] = mean(tmpS)
            sdValue[i,j] = sd(tmpS) / sqrt(length(tmpS))
        }
        
        cat(" & ")
        cat(mName[i], " ")
        for(j  in colT3){
            ss = sprintf("%0.3f",round(mValue[i,j],3))
            cat(" & ",ss )
        }
        cat(" \\\\\n ")
        ## row of standard error
        cat(" &  ") ## empty method name
        for(j  in colT3){
            ss = sprintf("(%0.3f)",round(sdValue[i,j],3))
            cat(" & ",ss )
        }
        cat(" \\\\\n")
        
    }
    return(mValue)
}


rowN = 5
colN = 7
colT2 = 1:7
colT3 = 1:7
specialList = c()

load("./data/simu1_1.RData")
res = get_mean_value(Simu1_1,  colT2, colT3)

# ----------- simu2 ------------
load("./data/simu2_1.RData")
res = get_mean_value(Simu2_1,  colT2, colT3)

# ----------- simu3 ------------


laGP_print = function(resMat){
    res1 = round(colMeans(resMat),3)
    res2 = round(apply(resMat, 2,sd) /sqrt(100),3)
    nn = length(res1)
    
    cat(" & laGP ")
    
    for(i in 1:nn){
        cat("& ",sprintf("%.3f",res1[i]),"")
    }
    
    cat("\\\\\n")
    cat(" &  ")
    for(i in 1:nn){
        cat("& (",sprintf("%.3f",res2[i]),")  ",sep = "")
    }
    cat("\\\\\n")
    
}


load("./data/simu3-1.RData")
res = get_mean_value(resultSF,  colT2, colT3)
load("./data/simu3_laGP.RData")
laGP_print(resMatCheap)

# ----------- simu4 ------------
load("./data/simu4-1.RData")
res = get_mean_value(resultSF2,  colT2, colT3)
load("./data/simu4_laGP.RData")
laGP_print(resMatCheap)


