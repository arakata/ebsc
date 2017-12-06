source("sourceHVK.R")

#simulate data
x=seq(0,1,length.out=256)
mu=(exp(1)^(10*x-5))/(1+exp(1)^(10*x-5))
#mu=doppler(x)+noise/10
noise=arima.sim(n=length(mu),list(ar=c(0.5)),sd=0.1)
y=mu+noise

#HVK
o<-HVKd(xvec=x,m1=NULL,m2=NULL,Y=y)
plot(y)
lines(o$yhat,lwd=2,col=4)

