source("sourceNason.R")

#simulate data
x=seq(0,1,length.out=256)
noise=arima.sim(n=length(x),list(ar=c(0.5)),sd=0.1)
#mu=doppler(x)+noise/10
mu=(exp(1)^(10*x-5))/(1+exp(1)^(10*x-5))
y=mu+noise

#Nason's Book
yhat<-nason(y)$yhat
plot(y)
lines(yhat,col=4,lwd=2)


