source("sourcelmm2.R")

#simulate data
x=seq(0,1,length.out=256)
#x=seq(0,1,length.out=200)
mu=(exp(1)^(10*x-5))/(1+exp(1)^(10*x-5))
#mu=doppler(x)
noise=arima.sim(n=length(x),list(ar=c(0.5)),sd=0.1)
y=mu+noise

#lmm2
yhat<-lmm2(y,c(1,0))$yhat
par(mfrow = c(1,1))
plot(x,y)
lines(x,yhat,col=4,lwd=2)