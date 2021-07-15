

ComputeThetaLoss = function(TrueSeqX, TrueSeqTheta, ThetaHat, ThetaSigma){
    deltaX = diff(TrueSeqX)[1]
    diff = TrueSeqTheta - ThetaHat
    L2Loss = sqrt(sum(diff^2) * deltaX)
    CIWidth = sum(ThetaSigma * 2 * 1.645) * deltaX
    
    
    UpperB = ThetaHat + 1.645 * ThetaSigma
    LowerB = ThetaHat - 1.645 * ThetaSigma
    coverB = LowerB < TrueSeqTheta & TrueSeqTheta < UpperB
    cm1_90 = sum(coverB) / length(coverB)
    width_90 = sum(UpperB - LowerB) * deltaX
    
    UpperB = ThetaHat + 1.959964 * ThetaSigma
    LowerB = ThetaHat - 1.959964 * ThetaSigma
    coverB = LowerB < TrueSeqTheta & TrueSeqTheta < UpperB
    cm1_95 = sum(coverB) / length(coverB)
    width_95 = sum(UpperB - LowerB) * deltaX
    
    UpperB = ThetaHat + 2.575829 * ThetaSigma
    LowerB = ThetaHat - 2.575829 * ThetaSigma
    coverB = LowerB < TrueSeqTheta & TrueSeqTheta < UpperB
    cm1_99 = sum(coverB) / length(coverB)
    width_99 = sum(UpperB - LowerB) * deltaX
    res = c(L2Loss, width_90,  cm1_90, width_95, cm1_95, width_99, cm1_99)
    return(res)
}


# Full curve coverage error, instead of pointwise
ComputeThetaLossWhole = function(TrueSeqX, TrueSeqTheta, ThetaHat, ThetaSigma){
    deltaX = diff(TrueSeqX)[1]
    diff = TrueSeqTheta - ThetaHat
    L2Loss = sqrt(sum(diff^2) * deltaX)
    CIWidth = sum(ThetaSigma * 2 * 1.645) * deltaX
    
    
    UpperB = ThetaHat + 1.645 * ThetaSigma
    LowerB = ThetaHat - 1.645 * ThetaSigma
    coverB = LowerB < TrueSeqTheta & TrueSeqTheta < UpperB
    cm1_90 = all(coverB)#sum(coverB) / length(coverB)
    width_90 = sum(UpperB - LowerB) * deltaX
    
    UpperB = ThetaHat + 1.959964 * ThetaSigma
    LowerB = ThetaHat - 1.959964 * ThetaSigma
    coverB = LowerB < TrueSeqTheta & TrueSeqTheta < UpperB
    cm1_95 = all(coverB)#sum(coverB) / length(coverB)
    width_95 = sum(UpperB - LowerB) * deltaX
    
    UpperB = ThetaHat + 2.575829 * ThetaSigma
    LowerB = ThetaHat - 2.575829 * ThetaSigma
    coverB = LowerB < TrueSeqTheta & TrueSeqTheta < UpperB
    cm1_99 = all(coverB)#sum(coverB) / length(coverB)
    width_99 = sum(UpperB - LowerB) * deltaX
    res = c(L2Loss, width_90,  cm1_90, width_95, cm1_95, width_99, cm1_99)
    return(res)
}
