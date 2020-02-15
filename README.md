# groupPLS
two step penalized group partial least squares

library(devtools) 
install_github("jialiwang1211/groupPLS")
library(groupPLS)

data<-read.csv("Spectra_eps.csv")
data<-data[complete.cases(data),]
groupPLS(data[,c(2,4:6)] ~ data[,7:dim(data)[2]], k=3)
