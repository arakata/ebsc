source("sourceWang.R")

#simulate data
x=seq(0,1,length.out=256)
mu=(exp(1)^(10*x-5))/(1+exp(1)^(10*x-5))
#mu=doppler(x)
plot(y,col=8)
for(i in 1:100){
print(paste("n=",i))
noise=arima.sim(n=length(x),list(ar=c(0.5)),sd=0.1)
y=mu+noise

#Cochrane-Orcut iterations (Wang)
o<-CochraneOrcut(num.knots=5,bdeg=2,y=y)
print(o$p.star)
#plot(y)
#lines(o$yhat,lwd=2,col=4)
lines(o$yhat,lwd=1,col=4)
}
