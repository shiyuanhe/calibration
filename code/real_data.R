rm(list = ls())
source("./code/core/loadAll.R")
source("./code/core/real_core.R")

mData = R.matlab::readMat("./data/data_new_Phy_Sim.mat")

real_x = mData$xa#[-c(1,3)]#[-c(1:2,14:17)]
real_yp_original = mData$PA#[-c(1,3)]#[-c(16,17)]#[-c(1:2,14:17)] 
#plot(real_x, real_yp_original)

emu_x = mData$PVA.a[,1]
emu_theta = mData$PVA.a[,2]
emu_y = mData$A.502[,2]
emuData = list(x = cbind(emu_x, emu_theta), y = emu_y)
emu_betaOpt = c(3.932682, -2.6915)

lowerB = 0.4
upperB = 1.2
nSample = length(real_x)
real_x = matrix(real_x, nSample)
TrueSeqX = real_x
predX = matrix(seq(lowerB, upperB, length.out = 100), nrow =100)
linkLowerB = 0.35
linkUpperB = 1.05



### Debias and plot
sel = real_x > 0.5 & real_x < 1.0
bias = mean(emu_y) - mean(real_yp_original[sel]) 
real_yp = real_yp_original + bias
library(ggplot2)

#pdf("./figure/real_data_debiased.pdf",height = 5, width = 6)
emu_data = data.frame(emu_x = emu_x, emu_y = emu_y, emu_theta = emu_theta)
real_data = data.frame(real_x = real_x, real_yp = real_yp, emu_theta = 1)
real_data2 = data.frame(real_x = real_x, real_yp = real_yp_original, emu_theta = 1)
ggplot(emu_data, aes(x = emu_x, y = emu_y, color = emu_theta)) + geom_point(shape = 17) +
    geom_point(aes(real_x, real_yp, fill = "debiased"), pch = 21,data = real_data,size = 2) +
    geom_point(aes(real_x, real_yp, fill = "original"), pch = 21,data = real_data2,size = 2) +
    scale_color_gradient(low="steelblue",high="white")+
    #geom_point(aes(real_x, yHat), data = cv_data, color = "red",alpha = 1) +
    #geom_point(aes(real_x, yHat), data = cv_data_const, color = "green",alpha = 1) +
    theme_bw() + xlab("X") + ylab("Y")+
    scale_fill_manual(values=c("red", "blue","blue"))
#dev.off()


#####  RKHS-Cubic #####
#leave-one-out result 
lambdaSeqCubic = (exp(seq(-4,4,length.out = 25)))
yHat = RKHS_OOB(c2 = 5, lambdaSeqCubic, getInitGamma_cal6)
(rkhsRes1 = computeErr(yHat))
# random leave-two-out result
yHat2 = RKHS_OOB_randOUT(c2 = 5, 100, lambdaSeqCubic, getInitGamma_cal6,oobSize = 2)
(rkhsRes2 = c(mean(yHat2/2), sd(yHat2/2)/sqrt(100)))
