# groupPLS
Perform two step penalized group partial least squares on spectroscopy data from a eucalypt forest. Spectral reflectance was measured for 69 healthy leaves, with spectral data between 350 - 2500 nm and 1 nm band spacing resolution. We predicted the four response variables: carotene, chlorophyll, EPS (xanthophyll epoxidation state) and weight using the full spectrum data.


library(devtools) \
install_github("jialiwang1211/groupPLS")\
library(groupPLS)

data<-read.csv("Spectra_eps.csv")\
data<-data[complete.cases(data),]\

groupPLS(data[,c(2,4:6)] ~ data[,7:dim(data)[2]], k=3)\
