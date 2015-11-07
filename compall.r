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
z=2
ttheta=0.4
Amax=0
a_lo=-wl/(1+i)+1.0e-2
a_hi=b+wl*(1-h)-Ph*h-1.0e-2

seq(a_lo,a_hi, 0.001)



Assetsolver(b,wl,wh,Ph,h,eeta,i,ppsi,ttheta,bbeta,5)

Assetsolver(0.03,wl,wh,Ph,h,eeta,i,ppsi,0.59,bbeta,5)





bgrid<-seq(1,5,(5/100))
tthetagrid<-seq(0.1,0.9,(1/100))
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


filled.contour(bgrid,tthetagrid,hopt, main="EEDU")


algo1=Assetsolver(bgrid[bb],wl,wh,Ph,0,eeta,i,ppsi,tthetagrid[tt],bbeta,0)
algo2=Assetsolver(bgrid[bb],wl,wh,Ph,1,eeta,i,ppsi,tthetagrid[tt],bbeta,0)
bbmin=1
bbstep=0.1
bbmax=5
tthetamin=0.1
tthetastep=0.01
tthetamax=0.9
aalphah=0.8

pphi=0.2
ggamma=0.4


##Parameter definitions
bbeta=0.971
ppsi=2
A=1
aalphah=0.5
z=1.5
pphi=2
ggamma=0.66
eeta=0.5
wl=1
wh=7
i=0.02
GENERAL=c(log(wl),log(wh),log(i))
Ph=1
bbmin=1
bbmax=5
bbstep=(5-1)/100
tthetamin=0.1
tthetamax=0.9
tthetastep=(0.9-0.1)/100
Amax=0
par=c(ppsi,bbeta,eeta,aalphah,z,pphi,ggamma,A)
start = proc.time()[3];
a=hres(par,GENERAL,Ph,Amax,tthetamin,tthetamax,tthetastep,
     bbmin,bbmax,bbstep)


finish = proc.time()[3] - start;
print(finish)




#Optimizing the likelihood function
finfun<- function(x){
  restotal=hres(par,x,Ph,Amax,tthetamin,tthetamax,tthetastep,
                bbmin,bbmax,bbstep)
  return(restotal)
} 



OPT=optim(GENERAL, finfun, method="BFGS", hessian=T, control=list(reltol=1.0e-16, maxit=500000))
paropt=OPT$par
