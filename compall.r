
#0. Housekeeping
#Cleanning everything

rm(list=ls())

#Importing the datasete
#install.packages("pkgKitten")
#install.packages("knitr")
#install.packages("rmarkdown")
#install.packages("Rcpp")
#install.packages("RcppArmadillo")
#install.packages("cppFunction")
#install.packages("BH")
#install.packages("RcppGSL")
library(RcppGSL)
library(RcppArmadillo)
library(Rcpp)
library(BH)

Rcpp::sourceCpp('Documents/Research/zarruk/Didris/compall.cpp')


bbeta=0.971
ppsi=2
eeta=0.5
wl=1
wh=7
i=0.02
Ph=1
h=1
b=2
ttheta=0.4

a_lo=-wl/(1+i)+1.0e-2
a_hi=b+wl*(1-h)-Ph*h-1.0e-2
F_root_ab(a_lo,b,wl,wh,Ph,h,eeta,i,ppsi,ttheta,bbeta)
F_root_ab(a_hi,b,wl,wh,Ph,h,eeta,i,ppsi,ttheta,bbeta)
seq(a_lo,a_hi, 0.001)
F_root_ab(-100000,b,wl,wh,Ph,h,eeta,i,ppsi,ttheta,bbeta)
F_root_ab(1000000,b,wl,wh,Ph,h,eeta,i,ppsi,ttheta,bbeta)
F_root_ab(10,b,wl,wh,Ph,h,eeta,i,ppsi,ttheta,bbeta)

F_root_ab(-0.8,b,wl,wh,Ph,h,eeta,i,ppsi,ttheta,bbeta)
F_root_ab(0.8,b,wl,wh,Ph,h,eeta,i,ppsi,ttheta,bbeta)


Assetsolver(b,wl,wh,Ph,h,eeta,i,ppsi,ttheta,bbeta,5)

