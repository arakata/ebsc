source("sourceNason.R")

#simulate data
x=seq(0,1,length.out=256)
mu=(exp(1)^(10*x-5))/(1+exp(1)^(10*x-5))
noise=arima.sim(n=length(x),list(ar=c(0.5)),sd=0.1)
y=mu+noise

#Nason's Book
yhat<-nason(y)$yhat
that<-nason(y)$t
plot(x,y)
lines(that,yhat,col=4,lwd=2)

# other example
n=10^3
tt <- sort(runif(n))
mu<-doppler(tt)
yy <-mu + arima.sim(n=n,list(ar=c(0.05)),sd=0.15)

yygrid <- makegrid(t=tt, y=yy)
yyirregwd <- irregwd(yygrid)
yy.thresh.sure <- threshold(yyirregwd, policy="sure", type="soft", dev=madmad) # policy sure is suppose to deal with correlation

yy.wr <- wr(yy.thresh.sure)
plot(tt, yy,col=8,type='l')
lines(tt,mu,col=4,lwd=3)
lines(yygrid$gridt, yy.wr, col=2,lwd=3)
