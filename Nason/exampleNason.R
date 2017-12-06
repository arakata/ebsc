source("sourceNason.R")

#simulate data
#x=seq(0,1,length.out=256)
#x=seq(0,1,length.out=200)
mu=(exp(1)^(10*x-5))/(1+exp(1)^(10*x-5))
#mu=doppler(x)
noise=arima.sim(n=length(x),list(ar=c(0.5)),sd=0.1)
y=mu+noise

#Nason's Book
yhat<-nason(y)$yhat
that<-nason(y)$t
plot(x,y)
lines(that,yhat,col=4,lwd=2)


