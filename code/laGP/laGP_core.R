library(laGP)
library(crs)

laGP_print = function(resMat){
    res1 = round(colMeans(resMat),3)
    res2 = round(apply(resMat, 2,sd) /sqrt(100),3)
    nn = length(res1)

    cat(" & laGP &")
    
    for(i in 1:nn){
        cat(sprintf("%.3f",res1[i]),"& ")
    }
    
    cat("\n")
    cat(" &  &")
    for(i in 1:nn){
        cat(" (",sprintf("%.3f",res2[i]),") & ",sep = "")
    }
    
}

# real Dataset
# emulation dataset
# The final prediction is made over fullX
# Output: the computed y at fullX and its sigma uncertainty
laGP_predict = function(real_x, real_yp, emuData, fullX){
    
    real_x = simuData$x
    real_yp  = simuData$y
    fullX = TrueSeqX
    
    sel1 = 1:xDim
    emu_x = emuData$x[,sel1,drop=FALSE]
    sel2 = (xDim+1):(xDim+ncalib)
    emu_theta = emuData$x[,sel2]
    emu_y = emuData$y
    XU = emuData$x
    da <- d <- darg(NULL, XU) # length scale parameter for 
    
    # overwrite 
    # d$start = 0.5 # real data 0.1
    # da = d
    
    #methods = FALSE to use the computer model M directly for fcalib and aGP
    if(!is.null(cheapModel)){
        methods = FALSE 
        M = cheapModel
    } 

    g <- garg(list(mle=TRUE), real_yp) # nugget parameter
    
    if(!is.null(hyperG)) g$start = hyperG
    if(!is.null(hyperD)) g$start = hyperD
    
    ## the calibration parameter, u
    if(ncalib == 2){
        theta_seq1 = runif(10,snomadr_lb[1], snomadr_ub[1])
        theta_seq2 = runif(10,snomadr_lb[2], snomadr_ub[2])
        uinit <- cbind(theta_seq1, theta_seq2)
    }
    
    if(ncalib==1) uinit = matrix(runif(10,snomadr_lb, snomadr_ub),ncol = 1)
    
    if(ncalib>2) stop("ncalib>2")
    
    llinit <- rep(NA, nrow(uinit))
    for(i in 1:nrow(uinit)) {
        # XU, Z emulator data
        # X, Y field (read) data (Y double size??)
        # M might be optional
        llinit[i] <- fcalib(uinit[i,], XU, emu_y, real_x, real_yp, 
                            da, d, g, beta.prior, 
                            methods, cheapModel, bias.est, nThread)
        # evalute the loglikelihood or posterior
    }
    
    
    
    ## now for a NOMAD search via snomadr in the crs package
    its <- 0
    o <- order(llinit)
    i <- 1
    out <- NULL
    while(its < 10) {
        cat("NOMAD start=", i, ", uinit=(", paste(uinit[o[i],], collapse=", "), ")\n", sep="")
        outi <- snomadr(fcalib, snomadr_n, snomadr_bbin, snomadr_bbout, 
                        x0=uinit[o[i],],
                        lb=snomadr_lb, ub=snomadr_ub, 
                        opts=opts, XU=XU, 
                        Z=emu_y, X=real_x, Y=real_yp, 
                        da=da, d=d, g=g, methods=methods, M=cheapModel, 
                        bias=bias.est, omp.threads=nThread, uprior=beta.prior, 
                        save.global=.GlobalEnv, verb=0)
        its <- its + outi$iterations
        ## keep the best
        if(is.null(out) || outi$objective < out$objective) out <- outi
        i <- i + 1;
    }
    
    u.hat <- outi$solution
    
    # g$start = 0.2
    # g$max = 2
    # d$start = 0.2
    if(!is.null(hyperG)) g$start = hyperG
    if(!is.null(hyperD)) g$start = hyperD
    ## now with the estimated u value
    nSample = nrow(real_x)
    X_uhat <- cbind(real_x, matrix(rep(u.hat, nSample), ncol=ncalib, byrow=TRUE))
    Mhat <- aGP.seq(XU, emu_y, X_uhat, da, methods, ncalib=ncalib,
                    omp.threads=nThread, verb=0, M = cheapModel)
    cmle <- discrep.est(real_x, real_yp, Mhat$mean, d, g, bias.est, FALSE)
    
    
    # Predict over the test set
    nTest = nrow(fullX)
    XSeq_uhat = cbind(fullX, matrix(rep(u.hat, nTest), ncol=ncalib, byrow=TRUE))
    Mhat.oos <- aGP.seq(XU, emu_y, XSeq_uhat, da, methods, ncalib=ncalib,
                        omp.threads=nThread, verb=0, M = cheapModel)
    YYm.pred <- predGP(cmle$gp, matrix(XSeq_uhat[,1], ncol = xDim))
    YY.pred <- YYm.pred$mean + Mhat.oos$mean
    ySigma <- sqrt(diag(YYm.pred$Sigma))
    # 
    # plot(TrueSeqX, TrueSeqY, type = "l")
    # lines(TrueSeqX, YY.pred, col = "red")
    # lines(TrueSeqX, YY.pred +ySigma, col = "red",lty = 2)
    # lines(TrueSeqX, YY.pred -ySigma, col = "red",lty = 2)
    
    return(cbind(YY.pred, ySigma))
}
