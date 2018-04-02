library(eBsc)

q.max = 6;
all_ns = c(100, 2^6, 2^7, 2^8)
Basis = list()
for (i_q in 1:q.max)
{
	Basis[[i_q]] = list()
    	for(i_n in 1:length(all_ns))
      	Basis[[i_q]][[i_n]] = drbasis(nn = all_ns[i_n], qq = i_q)
}
  
ns = c(2^7, 2^8) # chosen values of n to try out
n_ns = length(ns)
sigma = 0.1 # standard deviations for the noise

mu3 = Basis[[3]][[which( all_ns == max(all_ns))]]$eigenvectors %*% (c(rep(0,5), (pi*(3:(max(all_ns)-3)))^(-(3 + 0.1)))*(cos(2*(1:max(all_ns)))))
mu3 = mu3/(max(mu3)-min(mu3))

mu5 = Basis[[5]][[which( all_ns == max(all_ns))]]$eigenvectors %*% (c(rep(0,5), (pi*(3:(max(all_ns)-3)))^(-(5 + 0.1)))*(cos(2*(1:max(all_ns)))))
mu5 = mu5/(max(mu5)-min(mu5))

## functions to be considered
## make sure that functions' range has length 1 so that the results are on the same scale
function_names = c()
functions = list()
function_names = c(function_names, "a 3-smooth function")
functions[[1]] = function(x) predict(interpSpline(seq(0,1,length=ns[n_ns]), mu3, bSpline = T), x)$y
function_names = c(function_names, "a 5-smooth function")
functions[[2]] = function(x) predict(interpSpline(seq(0,1,length=ns[n_ns]), mu5, bSpline = T), x)$y
function_names = c(function_names, "a mixture of two Beta densities")
functions[[3]] = function(x) (6*dbeta(x,30,17)+4*dbeta(x,3,11))/10/3.35
function_names = c(function_names, "a mixture of three Beta densities")
functions[[4]] = function(x) (dbeta(x,10,5)+dbeta(x,7,7)+dbeta(x,5,10))/5.38
function_names = c(function_names, "a cosine function")
functions[[5]] = function(x) cos(2*pi*x+pi)/2
# function_names = c(function_names, "name for the function")
# functions[[4]] = function(x) ???
n_functions = length(functions)


## noise processes to be considered
ARMA_pars = list()
ARMA_pars[[1]] = list(ar = c(0.6), ma = c())
ARMA_pars[[2]] = list(ar = c(-0.5, 0.1), ma = c())
ARMA_pars[[3]] = list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488))
n_correlations = length(ARMA_pars)


## generate one dataset for each combination of function and correlation, and get estimates
n = 2^7
data = list()
truths = list()
result = list()
for(i_function in 1:n_functions)
{
  data[[i_function]] = list()
  truths[[i_function]] = functions[[i_function]](seq(0,1, length = n))
  for(i_correlation in 1:n_correlations)
  {
    data[[i_function]][[i_correlation]] = truths[[i_function]] + sigma * arima.sim(n = n, ARMA_pars[[i_correlation]], sd = 1)
    }
}


DATA.T<-cbind(truths[[1]],truths[[2]],truths[[3]],truths[[4]],truths[[5]])
write.csv2(DATA.T,file="truths.csv",row.names = FALSE)

#f1<-cbind(data[[1]][[1]],data[[1]][[2]],data[[1]][[3]])
#f2<-cbind(data[[2]][[1]],data[[2]][[2]],data[[2]][[3]])
#f3<-cbind(data[[3]][[1]],data[[3]][[2]],data[[3]][[3]])
#f4<-cbind(data[[4]][[1]],data[[4]][[2]],data[[4]][[3]])
#f5<-cbind(data[[5]][[1]],data[[5]][[2]],data[[5]][[3]])

#write.csv2(f1,file="simf1.csv",row.names = FALSE)
#write.csv2(f2,file="simf2.csv",row.names = FALSE)
#write.csv2(f3,file="simf3.csv",row.names = FALSE)
#write.csv2(f4,file="simf4.csv",row.names = FALSE)
#write.csv2(f5,file="simf5.csv",row.names = FALSE)

####
#f1
####

nf1c1<-NULL
for(i in 1:1000){
    nf1c1<-cbind(nf1c1,DATA.T[,1] + sigma *arima.sim(n = n, list(ar = c(0.6), ma = c()), sd = 1))
}

nf1c2<-NULL
for(i in 1:1000){
    nf1c2<-cbind(nf1c2,DATA.T[,1] + sigma *arima.sim(n = n, list(ar = c(-0.5, 0.1), ma = c()), sd = 1))
}

nf1c3<-NULL
for(i in 1:1000){
    nf1c3<-cbind(nf1c3,DATA.T[,1] + sigma *arima.sim(n = n, list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488)), sd = 1))
}

####
#f2
####

nf2c1<-NULL
for(i in 1:1000){
    nf2c1<-cbind(nf2c1,DATA.T[,2] + sigma *arima.sim(n = n, list(ar = c(0.6), ma = c()), sd = 1))
}

nf2c2<-NULL
for(i in 1:1000){
    nf2c2<-cbind(nf2c2,DATA.T[,2] + sigma *arima.sim(n = n, list(ar = c(-0.5, 0.1), ma = c()), sd = 1))
}

nf2c3<-NULL
for(i in 1:1000){
    nf2c3<-cbind(nf2c3,DATA.T[,2] + sigma *arima.sim(n = n, list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488)), sd = 1))
}


####
#f3
####

nf3c1<-NULL
for(i in 1:1000){
    nf3c1<-cbind(nf3c1,DATA.T[,3] + sigma *arima.sim(n = n, list(ar = c(0.6), ma = c()), sd = 1))
}

nf3c2<-NULL
for(i in 1:1000){
    nf3c2<-cbind(nf3c2,DATA.T[,3] + sigma *arima.sim(n = n, list(ar = c(-0.5, 0.1), ma = c()), sd = 1))
}

nf3c3<-NULL
for(i in 1:1000){
    nf3c3<-cbind(nf3c3,DATA.T[,3] + sigma *arima.sim(n = n, list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488)), sd = 1))
}


####
#f4
####

nf4c1<-NULL
for(i in 1:1000){
    nf4c1<-cbind(nf4c1,DATA.T[,4] + sigma *arima.sim(n = n, list(ar = c(0.6), ma = c()), sd = 1))
}

nf4c2<-NULL
for(i in 1:1000){
    nf4c2<-cbind(nf4c2,DATA.T[,4] + sigma *arima.sim(n = n, list(ar = c(-0.5, 0.1), ma = c()), sd = 1))
}

nf4c3<-NULL
for(i in 1:1000){
    nf4c3<-cbind(nf4c3,DATA.T[,4] + sigma *arima.sim(n = n, list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488)), sd = 1))
}


####
#f5
####

nf5c1<-NULL
for(i in 1:1000){
    nf5c1<-cbind(nf5c1,DATA.T[,5] + sigma *arima.sim(n = n, list(ar = c(0.6), ma = c()), sd = 1))
}

nf5c2<-NULL
for(i in 1:1000){
    nf5c2<-cbind(nf5c2,DATA.T[,5] + sigma *arima.sim(n = n, list(ar = c(-0.5, 0.1), ma = c()), sd = 1))
}

nf5c3<-NULL
for(i in 1:1000){
    nf5c3<-cbind(nf5c3,DATA.T[,5] + sigma *arima.sim(n = n, list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488)), sd = 1))
}

#######
#GRABAR
#######

write.csv2(nf1c1,file="simf1c1.csv",row.names = FALSE)
write.csv2(nf1c2,file="simf1c2.csv",row.names = FALSE)
write.csv2(nf1c3,file="simf1c3.csv",row.names = FALSE)

write.csv2(nf2c1,file="simf2c1.csv",row.names = FALSE)
write.csv2(nf2c2,file="simf2c2.csv",row.names = FALSE)
write.csv2(nf2c3,file="simf2c3.csv",row.names = FALSE)

write.csv2(nf3c1,file="simf3c1.csv",row.names = FALSE)
write.csv2(nf3c2,file="simf3c2.csv",row.names = FALSE)
write.csv2(nf3c3,file="simf3c3.csv",row.names = FALSE)

write.csv2(nf4c1,file="simf4c1.csv",row.names = FALSE)
write.csv2(nf4c2,file="simf4c2.csv",row.names = FALSE)
write.csv2(nf4c3,file="simf4c3.csv",row.names = FALSE)

write.csv2(nf5c1,file="simf5c1.csv",row.names = FALSE)
write.csv2(nf5c2,file="simf5c2.csv",row.names = FALSE)
write.csv2(nf5c3,file="simf5c3.csv",row.names = FALSE)

