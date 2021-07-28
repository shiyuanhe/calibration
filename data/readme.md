- `data_new_Phy_Sim.mat` is the Matlab format data for real data analysis.
- The other `RData` files are temporarily cached file to generate the `NumericalResults.pdf` file. 


```r
mData = R.matlab::readMat("./data/data_new_Phy_Sim.mat")

# Physical observation of size 17
mData$xa # input variable: PVA Amount
mData$PA # response variable: Modulus 

# Emulation data set of size 149
emu_x = mData$PVA.a[,1] # input variable: PVA Amount
emu_theta = mData$PVA.a[,2] # calibration variable: effectiveness
emu_y = mData$A.502[,2] # response variable: Modulus 
emuData = list(x = cbind(emu_x, emu_theta), y = emu_y)
```