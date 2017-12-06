####################################################################################################
## Numerics for correlated noise paper
####################################################################################################


## load libraries
library(eBsc)
library(mvtnorm)


## set working directories
try(setwd("/Users/pauloserra/Work/!Ongoing/Project T2/Numerics"), silent = TRUE) # Mac directory path
try(setwd(  "C:/Users/pserra/Work/!Ongoing/Project T2/Numerics"), silent = TRUE) # Windows directory path


## set parameters
set.seed(123)
cols = c("grey", "black", "red", "blue", "green", "orange", "purple", "ADD MORE") # list of colours so they are always used in the same order
lwds = 5 # thickness of lines in plots
q.max = 6 # maximal value for q


## options for parallel computing
reps = 10 # number of times that the simulations are repeated
WINDOWS = F # set to T if code is run on Windows machine

n_cores = detectCores() - 1 - (WINDOWS == T)
if(WINDOWS){
  cluster_type = "PSOCK"
  ## variables to export
  exports = c("methods", "method_1", "method_2", "method_3") # add more methods here
  exports = c(exports, "ns", "n_ns", "sigma", "functions", "n_functions", "ARMA_pars", "n_correlations", "simulation", "reps")
  exports = c(exports, "mu3", "mu5", "q.max", "n_methods", "approach_names")
  exports = c(exports, "interpSpline")
}else{
  cluster_type = "FORK"
}


## options for credible set simulation
alpha = 0.05 ## confidence level of the set is 1-alpha
reps_credible = 100 ## number of trajectories to sampel from the posterior


## options for saving / plotting
SAVE_PLOT = T # should plots be saved to disk? (if F then they are displayed)
RECOMPUTE = T # should all calculations be redone and saved to the disk? (if F then all computations are loaded from disk)




####################################################################################################
## functions implementing the different approaches

methods = list()
approach_names = NULL


compute_rho = function(n, ARMA_pars) ## REPLACE WITH SOMETHING THAT DOESN'T CHANGE THE ORDER OF THE EIGENVALUES
  return(eigen(toeplitz(ARMAacf(ar = ARMA_pars$ar, ma = ARMA_pars$ma, lag.max = n-1, pacf = FALSE)))$values)


approach_names = c(approach_names, "Our approach")
method_1 = function(data)
{# implementation of our estimator
  
  ## IMPLEMENTATION GOES HERE
  
  result = NULL # filled with dummy values for now
  result$data = data
  result$iterations = 2
  result$estimate  = data
  result$estimates = sapply(1:(result$iterations-1), function(x) data)
  result$lambda  = 1
  result$lambdas = rep(1, result$iterations-1)
  result$rho  = rep(1, length(data))
  result$rhos = sapply(1:(result$iterations-1), function(x) result$rho)
  result$R  = diag(length(data))
  result$S  = diag(length(data))
  result$sigmasq  = 1
  result$q  = 3
  result$qs = rep(3, result$iterations-1)
  result$Tq  = rep(0, q.max)
  result$Tqs = array(1, dim = c(q.max, result$iterations-1))
  
  return(result)
}
methods[[1]] = method_1


approach_names = c(approach_names, "Other approach")
method_2 = function(data)
{# implementation of the estimator of ????
  
  ## IMPLEMENTATION GOES HERE
  
  result = NULL # filled with dummy values for now
  result$estimate = data # WE MUST ALWAYS HAVE A FIELD $estimate
  
  return(result)
}
methods[[2]] = method_2


approach_names = c(approach_names, "Yet another approach")
method_3 = function(data)
{# implementation of the estimator of ????
  
  ## IMPLEMENTATION GOES HERE
  
  result = NULL # filled with dummy values for now
  result$estimate = data # WE MUST ALWAYS HAVE A FIELD $estimate
  
  return(result)
}
methods[[3]] = method_3


n_methods = length(methods)




####################################################################################################
## numerics for the real dataset


## load real data
data = read.csv("./data/bitcoin_price.csv") # just a placeholder so we have something nontrivial, we can pick another dataset
y = data$Close[rev(1:length(data$Close))]
n = length(y)
x = 1:n


## plot the real data
if(SAVE_PLOT) pdf(file = "./figures/dataset1_data_plot.pdf",width=50, height=30, family='Times', pointsize=36)
par(mfrow = c(1,1))
plot(x, y, type = "l", col = "blue", xlab = "Time", ylab = "US Dollars", lwd = lwds)
title("Closing price of Bitcoins in US Dollars")
if(SAVE_PLOT) dev.off()


## apply all methods
if(RECOMPUTE)
{
  result = list()
  for(i_method in 1:n_methods)
    result[[i_method]] = methods[[i_method]](y)
  save(result, file = "./results/estimates.csv")
}else
  load(file = "./results/estimates.csv")

## get estimates
estimates = array(NA, dim = c(n, n_methods))
for(i_method in 1:n_methods) estimates[,i_method] = result[[i_method]]$estimate


## plot results comparing data and all the smoothers
range_ = range(estimates)
if(SAVE_PLOT) pdf(file = "./figures/dataset1_results_estimators.pdf",width=50, height=30, family='Times', pointsize=36)
par(mfrow = c(1,1))
plot(x,y, ylim = range_, xlab = "Time", ylab = "US Dollars", col = "white", lwd = lwds)
title("Results of the estimation procedures")
points(x, y, col = cols[1])
for(i_method in 1:n_methods)
  lines(x, estimates[,i_method], col = cols[2+i_method], lwd = lwds)
legend("topleft", c("Observations", approach_names), 
       lwd = c(1, rep(lwds, n_methods)), col = cols[c(1,3:(3+n_methods-1))], 
       pch = c(1, rep(NA, n_methods)), lty = c(NA, rep(1,n_methods)))
if(SAVE_PLOT) dev.off()




##############################################################################################################################
## produce the credible set for real data

## get estimates of f, R, S, and the variance
f_hat = result[[1]]$estimate
R_hat = result[[1]]$R
S_hat = result[[1]]$S
s_hat = 30000*result[[1]]$sigmasq ## ??? throw away the multiplicative factor

## sample form the posterior
if(RECOMPUTE)
{
  fluctuations = rmvt(round(reps_credible/(1-alpha)), sigma = s_hat * ( S_hat %*% R_hat ), df = n+1, type = "shifted")
  fluctuations = fluctuations[order(apply(fluctuations^2, 1, sum)) <= reps_credible, ]
  save(fluctuations, file = "./results/fluctuations.csv")
}else
  load(file = "./results/fluctuations.csv")

## plot the credible set
if(SAVE_PLOT) pdf(file = "./figures/dataset1_results_credible.pdf",width=50, height=30, family='Times', pointsize=36)
par(mfrow = c(1,1))
plot(x, f_hat, type = "l", xlab = "Time", ylab = "US Dollars", col = "white", ylim = range(f_hat)+range(fluctuations))
title("Closing price of Bitcoins in US Dollars with credible set")
col_ = rgb(col2rgb(cols[2])[1]/255, col2rgb(cols[2])[2]/255, col2rgb(cols[2])[3]/255, 1.25/reps_credible)
for(i_line in 1:reps_credible)
  lines(x, f_hat + fluctuations[i_line, ], col = col_, lwd = lwds)
points(x, y, col = cols[1])
lines(x, f_hat, col = cols[4], lwd = lwds)
col_ = rgb(col2rgb(cols[2])[1]/255, col2rgb(cols[2])[2]/255, col2rgb(cols[2])[3]/255, 0.125)
legend("topleft", c("Observations","Posterior mean", "Credible set"), col = c(cols[c(1,4)], col_), 
       lwd = c(1, lwds, 3*lwds), pch = c(1, NA, NA), lty = c(NA, 1, 1))
if(SAVE_PLOT) dev.off()




####################################################################################################
## setup of numerics for simulated data

ns = c(50, 100, 200) # chosen values of n to try out
n_ns = length(ns)
sigma = 0.1 # standard deviations for the noise


## prepare the DR basis and DR-smooth functions
if(RECOMPUTE)
{
  Basis = list()
  for (i_q in 1:q.max)
  {
    Basis[[i_q]] = list()
    for(i_n in 1:n_ns)
      Basis[[i_q]][[i_n]] = drbasis(nn = ns[i_n], qq = i_q)
  }
  save(Basis, file = "./results/basis.csv")
}else
  load(file = "./results/basis.csv")


mu3 = Basis[[3]][[n_ns]]$eigenvectors %*% 
        (c(rep(0,5), (pi*(3:(ns[n_ns]-3)))^(-(3 + 0.1)))*(cos(2*(1:ns[n_ns]))))
mu3 = mu3/(max(mu3)-min(mu3))

mu5 = Basis[[5]][[n_ns]]$eigenvectors %*% 
        (c(rep(0,5), (pi*(3:(ns[n_ns]-3)))^(-(5 + 0.1)))*(cos(2*(1:ns[n_ns]))))
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


## plot the chosen functions
grid = seq(0,1, length = 1000)
if(SAVE_PLOT) pdf(file = "./figures/simulations_functions.pdf",width=30, height=50, family='Times', pointsize=36)
par(mfrow = c(n_functions,1))
for(i_function in 1:n_functions)
  plot(grid, functions[[i_function]](grid), col = cols[2+i_function], lwd = lwds, type = "l",
       main = paste0("Function ", i_function, ": ", function_names[i_function]), xlab = "x", 
       ylab = expression(f[i](x))) ## ???
if(SAVE_PLOT) dev.off()


## noise processes to be considered
ARMA_pars = list()
ARMA_pars[[1]] = list(ar = c(0.6), ma = c())
ARMA_pars[[2]] = list(ar = c(-0.5, 0.1), ma = c())
ARMA_pars[[3]] = list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488))
n_correlations = length(ARMA_pars)


## plot ACF from the chosen correlations
n = 50
grid = 0:(n-1)
if(SAVE_PLOT) pdf(file = "./figures/simulations_correlations.pdf",width=30, height=50, family='Times', pointsize=36)
par(mfrow = c(n_correlations,1))
for(i_correlation in 1:n_correlations)
  plot(grid, ARMAacf(ar = ARMA_pars[[i_correlation]]$ar, ma = ARMA_pars[[i_correlation]]$ma, 
                     lag.max = n-1, pacf = FALSE), type = "l", lwd = lwds, col = cols[2],
                      main = paste0("Autocorrelation function ", i_correlation),
                      xlab = "Lag", ylab = "Autocorrelation")
if(SAVE_PLOT) dev.off()










##############################################################################################################################
## produce the credible set for simulated data


## generate one dataset for each combination of function and correlation
n = 50
x = seq(0,1, length = n)
data = list()
truths = list()
for(i_function in 1:n_functions)
{
  data[[i_function]] = list()
  truths[[i_function]] = functions[[i_function]](seq(0,1, length = n))
  for(i_correlation in 1:n_correlations)
  {
    data[[i_function]][[i_correlation]] = truths[[i_function]] + sigma * arima.sim(n = n, ARMA_pars[[i_correlation]], sd = 1)
  }
}


## get estimates of f, R, S, and the variance for each dataset
result = list()
for(i_function in 1:n_functions)
{
  result[[i_function]] = list()
  for(i_correlation in 1:n_correlations)
  {
    temp = methods[[1]](data[[i_function]][[i_correlation]])
    result[[i_function]][[i_correlation]] = 1
    result[[i_function]][[i_correlation]]$f_hat = temp$estimate
    result[[i_function]][[i_correlation]]$R_hat = temp$R
    result[[i_function]][[i_correlation]]$S_hat = temp$S
    result[[i_function]][[i_correlation]]$s_hat = 3*sigma^2*temp$sigmasq ## ??? throw away the multiplicative factor
  }
}## fix this ???


## sample form the posterior
if(RECOMPUTE)
{
  fluctuations = list()
  for(i_function in 1:n_functions)
  {
    fluctuations[[i_function]] = list()
    for(i_correlation in 1:n_correlations)
    {
      R_hat = result[[i_function]][[i_correlation]]$R_hat
      S_hat = result[[i_function]][[i_correlation]]$S_hat
      s_hat = result[[i_function]][[i_correlation]]$s_hat
      
      temp = rmvt(round(reps_credible/(1-alpha)), sigma = s_hat * ( S_hat %*% R_hat ), df = n+1, type = "shifted")
      fluctuations[[i_function]][[i_correlation]] = temp[order(apply(temp^2, 1, sum)) <= reps_credible, ]
    }
  }
  save(fluctuations, file = "./results/fluctuations_simulations.csv")
}else
  load(file = "./results/fluctuations_simulations.csv")


## plot the credible set
if(SAVE_PLOT) pdf(file = "./figures/simulations_results_credible.pdf",width=50, height=30, family='Times', pointsize=36)
par(mfrow = c(n_functions,n_correlations))
for(i_function in 1:n_functions)for(i_correlation in 1:n_correlations)
{
  f_hat = result[[i_function]][[i_correlation]]$f_hat
  plot(x, f_hat, type = "l", xlab = "Time", ylab = "US Dollars", col = "white", 
       ylim = range(f_hat)+range(fluctuations[[i_function]][[i_correlation]]))
  title(paste0("Function ", i_function, " with correlation structure ", i_correlation))
  col_ = rgb(col2rgb(cols[2])[1]/255, col2rgb(cols[2])[2]/255, col2rgb(cols[2])[3]/255, 1.25/reps_credible)
  for(i_line in 1:reps_credible)
    lines(x, f_hat + fluctuations[[i_function]][[i_correlation]][i_line, ], col = col_, lwd = lwds)
  points(x, data[[i_function]][[i_correlation]], col = cols[1])
  lines(x, f_hat, col = cols[2+i_function], lwd = lwds)
  lines(x, truths[[i_function]], col = cols[2], lwd = lwds/2)
  col_ = rgb(col2rgb(cols[2])[1]/255, col2rgb(cols[2])[2]/255, col2rgb(cols[2])[3]/255, 0.125)
  legend("bottomleft", c("Observations", "True function", "Posterior mean", "Credible set"), col = c(cols[c(1,2,2+i_function)], col_), 
         lwd = c(1, lwds/2, lwds, 3*lwds), pch = c(1, NA, NA, NA), lty = c(NA,1, 1, 1))
}
if(SAVE_PLOT) dev.off()




####################################################################################################
## run numerics for simulated data


## function that applies all methods to all combinations of sample size, regression function, and correlation structure
simulation = function(dummy)
{
  results = list()
  for(i_n in 1:n_ns)
  {
    results[[i_n]] = list()
    for(i_function in 1:n_functions)
    {
      results[[i_n]][[i_function]] = list()
      for(i_correlation in 1:n_correlations)
      {
        y = functions[[i_function]](seq(0,1, length = ns[i_n])) + sigma * arima.sim(n = ns[i_n], ARMA_pars[[i_correlation]], sd = 1) # generate data
        
        result = list()
        for(i_method in 1:n_methods) # apply all methods to estimate the regression function
          result[[i_method]] = methods[[i_method]](y)
        
        results[[i_n]][[i_function]][[i_correlation]] = result
      }
    }
  }
  return(results)
}


## run the simulations
if(RECOMPUTE)
{
  st = Sys.time()
  try(stopCluster(cl), silent = TRUE)
  cl = makeCluster(n_cores, type = cluster_type)
  if(WINDOWS) clusterExport(cl, exports)
  results = parLapply(cl,1:reps, simulation)
  stopCluster(cl)
  tt = Sys.time() - st
  print(tt)
  save(results, file = "./results/results.csv")
}else
  load(file = "./results/results.csv")




####################################################################################################
## summarise our results


## plot the number of iterations
if(RECOMPUTE)
{
  extracted = array(NA, dim = c(n_functions, n_correlations, n_ns, reps))
  for(i_n in 1:n_ns)for(i_function in 1:n_functions)for(i_correlation in 1:n_correlations)
    extracted[i_function, i_correlation, i_n, ] = sapply(1:reps, function(i_rep) 
      results[[i_rep]][[i_n]][[i_function]][[i_correlation]][[1]]$iterations)
  save(extracted, file = "./results/results_iterations.csv")
}else
  load(file = "./results/results_iterations.csv")


if(SAVE_PLOT) pdf(file = "./figures/simulations_iterations.pdf",width=50, height=30, family='Times', pointsize=36)
par(mfrow = c(n_functions, n_correlations))
for(i_function in 1:n_functions)for(i_correlation in 1:n_correlations)
{
  title_text = paste0("Function ", i_function, " with correlation structure ", i_correlation)
  boxplot(names = ns, t(extracted[i_function, i_correlation, , ]), main = title_text, xlab = "n", ylab = "Iterations")
}
if(SAVE_PLOT) dev.off()


## plot the logarithms of the first lambdas and of the lambdas after the last iteration
if(RECOMPUTE)
{
  extracted = array(NA, dim = c(n_functions, n_correlations, n_ns, 2 * reps))
  for(i_n in 1:n_ns)for(i_function in 1:n_functions)for(i_correlation in 1:n_correlations)
  {
    extracted[i_function, i_correlation, i_n, ] = c(
      sapply(1:reps, function(i_rep) results[[i_rep]][[i_n]][[i_function]][[i_correlation]][[1]]$lambdas[1]), 
      sapply(1:reps, function(i_rep) results[[i_rep]][[i_n]][[i_function]][[i_correlation]][[1]]$lambda))
  }
  save(extracted, file = "./results/results_lambdas.csv")
}else
    load(file = "./results/results_lambdas.csv")
extracted = log(extracted) # to plot the logarithms of the lambdas instead of the lambdas themselves

## ADJUST THE RANGE OF THE PLOT SO THEY ARE ALL THE SAME
if(SAVE_PLOT) pdf(file = "./figures/simulations_lambdas.pdf",width=50, height=30, family='Times', pointsize=36)
par(mfrow = c(n_functions, n_correlations))
for(i_function in 1:n_functions)for(i_correlation in 1:n_correlations)
{
  df = data.frame(id = c(rep("first", reps), rep("last", reps)), t(extracted[i_function, i_correlation, , ]))
  
  title_text = paste0("Function ", i_function, " with correlation structure ", i_correlation)
  boxplot(df[,-1], xlim = c(0.5, ncol(df[,-1])+0.5), boxfill = "white", border = "white", names = ns, main = title_text,
          xlab = "n", ylab = "Logarithm of smoothing parameter")
  boxplot(df[which(df$id=="first"),-1], xaxt = "n", add = TRUE, boxfill= "grey",  boxwex=0.25, at = 1:ncol(df[,-1]) - 0.15)    
  boxplot(df[which(df$id=="last"), -1], xaxt = "n", add = TRUE, boxfill= "black", boxwex=0.25, at = 1:ncol(df[,-1]) + 0.15)
  if(i_function == 1 & i_correlation == 1)
    legend("topright", c("First iteration", "Last iteration"), col = c("grey", "black"), lty = 1, lwd = lwds)
}
if(SAVE_PLOT) dev.off()


## plot the MSEs of the estimates of the regression function
if(RECOMPUTE)
{
  extracted = array(NA, dim = c(n_functions, n_correlations, n_ns, 2 * reps))
  for(i_n in 1:n_ns)for(i_function in 1:n_functions)for(i_correlation in 1:n_correlations)
  {
    extracted[i_function, i_correlation, i_n, ] = c(
      sapply(1:reps, function(i_rep)
        mean( (results[[i_rep]][[i_n]][[i_function]][[i_correlation]][[1]]$estimates[,1] - functions[[i_function]](seq(0, 1, length = ns[i_n])))^2 ) ), 
      sapply(1:reps, function(i_rep)
        mean( (results[[i_rep]][[i_n]][[i_function]][[i_correlation]][[1]]$estimate      - functions[[i_function]](seq(0, 1, length = ns[i_n])))^2 ) ) )
  }
  save(extracted, file = "./results/results_MSEs_ours.csv")
}else
  load(file = "./results/results_MSEs_ours.csv")


range_ = range(extracted)
if(SAVE_PLOT) pdf(file = "./figures/simulations_MSEs.pdf",width=50, height=30, family='Times', pointsize=36)
par(mfrow = c(n_functions, n_correlations))
for(i_function in 1:n_functions)for(i_correlation in 1:n_correlations)
{
  df = data.frame(id = c(rep("first", reps), rep("last", reps)), t(extracted[i_function, i_correlation, , ]))
  
  title_text = paste0("Function ", i_function, " with correlation structure ", i_correlation)
  boxplot(df[,-1], xlim = c(0.5, ncol(df[,-1])+0.5), boxfill = "white", border = "white", names = ns, main = title_text, ylim = range_,
          xlab = "n", ylab = "Mean squared error")
  boxplot(df[which(df$id=="first"),-1], xaxt = "n", add = TRUE, border= "grey",  boxwex=0.25, at = 1:ncol(df[,-1]) - 0.15)    
  boxplot(df[which(df$id=="last"), -1], xaxt = "n", add = TRUE, border= "black", boxwex=0.25, at = 1:ncol(df[,-1]) + 0.15)
}
if(SAVE_PLOT) dev.off()


## plot the AEs of the estimates of rho
if(RECOMPUTE)
{
  extracted = array(NA, dim = c(n_functions, n_correlations, n_ns, 2 * reps))
  for(i_n in 1:n_ns)for(i_correlation in 1:n_correlations)
  {
    temp = compute_rho(ns[i_n], ARMA_pars[[i_correlation]])
    for(i_function in 1:n_functions)
    {
      extracted[i_function, i_correlation, i_n, ] = c(
        sapply(1:reps, function(i_rep)
          max( abs(results[[i_rep]][[i_n]][[i_function]][[i_correlation]][[1]]$rhos[,1] - temp) ) ), 
        sapply(1:reps, function(i_rep)
          max( abs(results[[i_rep]][[i_n]][[i_function]][[i_correlation]][[1]]$rho      - temp) ) ) )
    }
  }
  save(extracted, file = "./results/results_AEs.csv")
}else
    load(file = "./results/results_AEs.csv")


range_ = range(extracted)
if(SAVE_PLOT) pdf(file = "./figures/simulations_AEs.pdf",width=50, height=30, family='Times', pointsize=36)
par(mfrow = c(n_functions, n_correlations))
for(i_function in 1:n_functions)for(i_correlation in 1:n_correlations)
{
  df = data.frame(id = c(rep("first", reps), rep("last", reps)), t(extracted[i_function, i_correlation, , ]))
  
  title_text = paste0("Function ", i_function, " with correlation structure ", i_correlation)
  boxplot(df[,-1], xlim = c(0.5, ncol(df[,-1])+0.5), boxfill = "white", border = "white", names = ns, main = title_text, ylim = range_,
          xlab = "n", ylab = "Absolute error")
  boxplot(df[which(df$id=="first"),-1], xaxt = "n", add = TRUE, border= "grey",  boxwex=0.25, at = 1:ncol(df[,-1]) - 0.15)    
  boxplot(df[which(df$id=="last"), -1], xaxt = "n", add = TRUE, border= "black", boxwex=0.25, at = 1:ncol(df[,-1]) + 0.15)
}
if(SAVE_PLOT) dev.off()


## plot the values of q that we obtain
if(RECOMPUTE)
{
  extracted = array(NA, dim = c(n_functions, n_correlations, n_ns, 2 * reps))
  for(i_n in 1:n_ns)for(i_function in 1:n_functions)for(i_correlation in 1:n_correlations)
  {
    extracted[i_function, i_correlation, i_n, ] = c(
      sapply(1:reps, function(i_rep) results[[i_rep]][[i_n]][[i_function]][[i_correlation]][[1]]$qs[1]), 
      sapply(1:reps, function(i_rep) results[[i_rep]][[i_n]][[i_function]][[i_correlation]][[1]]$q))
  }
  save(extracted, file = "./results/results_qs.csv")
}else
  load(file = "./results/results_qs.csv")


## ADJUST THE RANGE OF THE PLOT SO THEY ARE ALL THE SAME
if(SAVE_PLOT) pdf(file = "./figures/simulations_qs.pdf",width=50, height=30, family='Times', pointsize=36)
par(mfrow = c(n_functions, n_correlations))
for(i_function in 1:n_functions)for(i_correlation in 1:n_correlations)
{
  df = data.frame(id = c(rep("first", reps), rep("last", reps)), t(extracted[i_function, i_correlation, , ]))
  
  title_text = paste0("Function ", i_function, " with correlation structure ", i_correlation)
  boxplot(df[,-1], xlim = c(0.5, ncol(df[,-1])+0.5), boxfill = "white", border = "white", names = ns, main = title_text,
          xlab = "n", ylab = "q")
  boxplot(df[which(df$id=="first"),-1], xaxt = "n", add = TRUE, boxfill= "grey",  boxwex=0.25, at = 1:ncol(df[,-1]) - 0.15)    
  boxplot(df[which(df$id=="last"), -1], xaxt = "n", add = TRUE, boxfill= "black", boxwex=0.25, at = 1:ncol(df[,-1]) + 0.15)
}
if(SAVE_PLOT) dev.off()




####################################################################################################
## compare our results with the competition


## Compute the MSEs of every combination of function, correlation, ns, and methods
if(RECOMPUTE)
{
  extracted = array(NA, dim = c(n_functions, n_correlations, n_ns, reps, n_methods))
  for(i_n in 1:n_ns)for(i_function in 1:n_functions)for(i_correlation in 1:n_correlations)for(i_method in 1:n_methods)
  {
    extracted[i_function, i_correlation, i_n, , i_method] =
      sapply(1:reps, function(i_rep)
        mean( (results[[i_rep]][[i_n]][[i_function]][[i_correlation]][[i_method]]$estimate - functions[[i_function]](seq(0, 1, length = ns[i_n])))^2 ) )
  }
  save(extracted, file = "./results/results_MSEs_all.csv")
}else
    load(file = "./results/results_MSEs_all.csv")


## here the MSEs are stacked to plot the MSEs in boxplots
if(RECOMPUTE)
  {
  extracted_stacked = array(NA, dim = c(n_functions, n_correlations, n_ns, n_methods * reps))
  for(i_n in 1:n_ns)for(i_function in 1:n_functions)for(i_correlation in 1:n_correlations)for(i_method in 1:n_methods)
  {
    extracted_stacked[i_function, i_correlation, i_n, 1:reps + (i_method-1)*reps] =
      sapply(1:reps, function(i_rep)
        mean( (results[[i_rep]][[i_n]][[i_function]][[i_correlation]][[i_method]]$estimate - functions[[i_function]](seq(0, 1, length = ns[i_n])))^2 ) )
  }
  save(extracted_stacked, file = "./results/results_MSEs_all_stacked.csv")
}else
  load(file = "./results/results_MSEs_all_stacked.csv")
  

## plot the MSEs in boxplots
range_ = range(extracted_stacked)
if(SAVE_PLOT) pdf(file = "./figures/simulations_comparison.pdf",width=50, height=30, family='Times', pointsize=36)
par(mfrow = c(n_functions, n_correlations))
for(i_function in 1:n_functions)for(i_correlation in 1:n_correlations)
{
  labs = sort(rep(1:n_methods,reps))
  df = data.frame(id = labs, t(extracted_stacked[i_function, i_correlation, , ]))
  
  title_text = paste0("Function ", i_function, " with correlation structure ", i_correlation)
  boxplot(df[,-1], xlim = c(0.5, ncol(df[,-1])+0.5), boxfill = "white", border = "white", names = ns, main = title_text, ylim = range_, xlab = "n", ylab = "Mean squared error")
  for(i_method in 1:n_methods)
    boxplot(df[which(df$id==i_method),-1], xaxt = "n", add = TRUE, border= cols[2+i_method],  boxwex=0.5/n_methods, at = 1:ncol(df[,-1]) - 0.15 + 0.6/n_methods * (i_method-1))
  if(i_function == 1 & i_correlation == 1)
    legend("topright", approach_names, col = cols[2+1:n_methods], lty = 1, lwd = lwds)
}
if(SAVE_PLOT) dev.off()




####################################################################################################
## END OF RUN ######################################################################################
####################################################################################################