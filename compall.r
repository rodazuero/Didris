#In this file I will try to fix the behavioral model
#The issue I have is with the aalpha41m and also with the
#variance of the behavior-> The closer to zero the better fit
#we get for some reason. In Behavioral4.cpp we checked all the functions

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
Amax=0
a_lo=-wl/(1+i)+1.0e-2
a_hi=b+wl*(1-h)-Ph*h-1.0e-2

seq(a_lo,a_hi, 0.001)



Assetsolver(b,wl,wh,Ph,h,eeta,i,ppsi,ttheta,bbeta,5)

Assetsolver(0.03,wl,wh,Ph,h,eeta,i,ppsi,0.59,bbeta,5)





bgrid<-seq(1,5,(5/50))
tthetagrid<-seq(0.1,0.9,(1/50))
SIZE=length(bgrid)
mat<-matrix(list(), nrow=length(bgrid), ncol=length(tthetagrid))
mat0<-matrix(list(), nrow=length(bgrid), ncol=length(tthetagrid))
hopt<-matrix(list(), nrow=length(bgrid), ncol=length(tthetagrid))
par=c(bbeta,eeta,ppsi)
for (bb in  1:length(bgrid)){
  for (tt in 1:length(tthetagrid)){
    mat0[bb,tt]=Assetsolver(bgrid[bb],wl,wh,Ph,0,eeta,i,ppsi,tthetagrid[tt],bbeta,0)
    mat[bb,tt]=Assetsolver(bgrid[bb],wl,wh,Ph,1,eeta,i,ppsi,tthetagrid[tt],bbeta,0)
    hopt[bb,tt]=Hopt(par,bgrid[bb],tthetagrid[tt],wl,wh,i,Ph,Amax)
  }
}

F_root_ab(-0.54234,wl,0,Ph,ppsi,bbeta,0.1,wh,eeta,1,i)

filled.contour(bgrid,tthetagrid,hopt, main="EEDU")


algo1=Assetsolver(bgrid[bb],wl,wh,Ph,0,eeta,i,ppsi,tthetagrid[tt],bbeta,0)
algo2=Assetsolver(bgrid[bb],wl,wh,Ph,1,eeta,i,ppsi,tthetagrid[tt],bbeta,0)


