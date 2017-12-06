#
# Generate some values in (0,1), then sort them (for plotting)
#
tt <- sort(runif(100))
#
# Now evaluate the \code{\link{doppler}} function and add
# some noise.
#
yy <- doppler(tt) + rnorm(100, 0, 0.15)
#
# Now make the grid with this data
#
yygrid <- makegrid(t=tt, y=yy)
#
# Jolly good. Now let's take the wavelet transform of this gridded data.
# Note that we have to use the \code{\link{irregwd}} function
# of the gridded data as it computes the variances of the coefficients
# as well as the coefficients themselves.
#
yyirregwd <- irregwd(yygrid,family="DaubExPhase")

#
# You might want to plot the coefficients
#
# If you want to see the actual coefficients you have to first convert
# the class of the yyirregwd object to a wd object and then use
# \code{\link{plot.wd}} like this
#
yyirregwd2 <- yyirregwd
class(yyirregwd2) <- "wd"
## Not run: plot(yyirregwd2)
#
# If you want to see the variance factors (essentially the coefficient
# variances divided by the overall variance). Then just use
# \code{\link{plot.irregwd}}
#
## Not run: plot(yyirregwd)
#
# Ok. So you've seen the coefficients. Now let's do some thresholding.
#
yy.thresh.sure <- threshold(yyirregwd, policy="sure", type="soft", dev=madmad)
#
# And now do the reconstruct
#
yy.wr <- wr(yy.thresh.sure)
#
# And you can even plot the answer on the new grid!
#
## Not run: 
plot(yygrid$gridt, yy.wr, type="l")
#
# And superimpose the original data!
#
## Not run: points(tt, yy)
points(tt, yy)
#
# This is sort of \code{Doppler} like!

#WaveThresh

#######################################
#According to the wavethresh documentation, policy=sure, is an alternative to deal
#with correllated noise. Donoho and Johnstone, 1995 or Johnstone and Silverman, 1997

n=10^2
tt <- sort(runif(n))
#tt <- seq(0.001953125,0.9980469,length.out=n)
mu<-doppler(tt)
yy <-mu + arima.sim(n=100,list(ar=c(0.5)),sd=0.25)
#yy <-mu + rnorm(n, 0, 0.35)

yygrid <- makegrid(t=tt, y=yy)
yyirregwd <- irregwd(yygrid)
yy.thresh.sure <- threshold(yyirregwd, policy="sure", type="soft", dev=madmad) 

yy.wr <- wr(yy.thresh.sure)
plot(yygrid$gridt, yy.wr, type="l")
points(tt, yy)
lines(tt,mu,col=4,lwd=3)

#####################
v <- DJ.EX()
x <- (1:1024)/1024 # Define X coordinates too
plot(x, v$bumps, type="l", ylab="Bumps")
ssig <- sd(v$bumps) # Bumps sd
SNR <- 2 # Fix our SNR
# Work out what the variance of the noise is...
sigma <- ssig/SNR
e <- rnorm(1024, mean=0, sd=sigma)
y <- v$bumps + e
plot(x, y, type="l", ylab="Noisy bumps")
#
# Plot wd of bumps
#
xlv <- seq(from=0, to=1.0, by=0.2)
bumpswd <- wd(v$bumps)
plot(bumpswd, main="", sub="",xlabvals=xlv*512, xlabchars=as.character(xlv),xlab="x")
#
# Plot wd of noisy bumps for comparison
#
ywd <- wd(y)
plot(ywd, main="", sub="",xlabvals=xlv*512, xlabchars=as.character(xlv),xlab="x")
FineCoefs <- accessD(ywd, lev=nlevels(ywd)-1)
sigma <- mad(FineCoefs)
utDJ <- sigma*sqrt(2*log(1024))
ywdT <- threshold(ywd, policy="manual", value=utDJ)
ywr <- wr(ywdT)
plot(x, ywr, type="l")
lines(x, v$bumps, lty=2)
ywdcvT <- threshold(ywd, policy="cv", dev=madmad)
ywrcv <- wr(ywdcvT)
plot(x, ywrcv, type="l", xlab="x", ylab="Cross-val.Estimate")
lines(x, v$bumps, lty=2)

#
# Load the EbayesThresh library
#
library("EbayesThresh")
#
# Threshold the noisy wavelet coefficients using EbayesThresh
#
ywdEBT <- ebayesthresh.wavelet(ywd)
#
# Do the inverse transform on the shrunk coefficients
#
ywrEB <- wr(ywdEBT)
#
# Plot the reconstruction over the original
#
x <- (1:1024)/1024
plot(x, ywrEB, type="l", xlab="x",ylab="EBayesThresh Estimate")
lines(x, v$bumps, lty=2)

#
# Compute NDWT (packet-ordered)
#
ywst <- wst(y)
#
# Plot it
#
plot(ywst, main="", sub="")
#
# Access the fine scale coefficients and compute
# universal threshold
#
FineWSTCoefs <- accessD(ywst, lev=nlevels(ywd)-1)
sigmaWST <- mad(FineWSTCoefs)
utWSTDJ <- sigmaWST*sqrt(2*log(1024))
#
# Threshold (default number of levels) using the
# universal threshold
#
ywstT <- threshold(ywst, policy="manual", value=utWSTDJ)
plot(ywstT, scaling="by.level", main="", sub="")

yABuv <- AvBasis(ywstT)
yl <- range(c(yABuv, v$bumps))
plot(x, yABuv, type="l", xlab="x",ylab="TI-universal Estimate", ylim=yl)
lines(x, v$bumps, lty=2)

#
# Create space for recording performance for
# each shifted basis. There is one shift for
# each element of y
#
rss <- rep(0, length(y))
#
# For each shift, i, first compute the node
# vector for that shift (which defines the
# basis). Then invert ywstT using the
# packets defined by the node vector.
# Form and store the measure of performance
#
for(i in 1:length(y)){
thenv <- numtonv(i-1, nlevels(ywstT))
therecon <- InvBasis(ywstT, nv=thenv)
rss[i] <- sqrt(sum( (therecon - bumps)^2))
}
# Report performance for the standard
# wavelet basis, and the best one.
#
cat("Standard wavelet basis RSS is ", rss[1], "\n")
cat("Best wavelet basis RSS is ", min(rss), "\n")
#
# Plot the performances
#
plot(1:length(y), rss, type="l",xlab="Basis Number",ylab="True Sum of Squares Error")

cmws <- cthresh(y, details=TRUE)
plot(cmws$data.wd)
plot(cmws$thr.wd)

yl <- range(c(Re(cmws$estimate), v$bumps))
plot(x, Re(cmws$estimate), type="l", xlab="x",ylab="Complex MW Estimate (Real)", ylim=yl)
lines(x, v$bumps, lty=2)

x <- c(rep(0,8), rep(1,16), rep(0,8))
xwd <- wd(x, filter.number=4)
round(accessD(xwd, level=4),3)

#
# Generate ARMA noise
#
eps <- arima.sim(n=1024, model=list(ar=0.99, ma=1))
#
# Create scaled noise with corrected SNR and
# then noisy signal
#
sigma<-0.1
eps <- sigma*eps/sqrt(var(eps))
y <- v$bumps + eps

plot(x, y, type = "l", ylab = "Correlated Noise bumps")
#
# Take wavelet transform #
ywd <- wd(y)
#
# Threshold BY LEVEL, but return the threshold
# value so we can print it out #
ywdT <- threshold(ywd, by.level=TRUE, policy="universal",return.thresh=TRUE)
print(ywdT)
#
# Now actually apply the threshold and invert
#
ywr <- wr(threshold(ywd, by.level=TRUE, policy="universal")) #
# Plot the denoised version and the original
#
yl <- range(c(ywr, v$bumps))
plot(x, ywr, ylim=yl, type="l")
lines(x, v$bumps, lty=2)

#############################################################
#BOOK by NASON: Wavelet Methods in Statistics with R page 135
#############################################################

v <- DJ.EX()
x <- (1:1024)/1024 # Define X coordinates too
plot(x, v$bumps, type="l", ylab="Bumps")
ssig <- sd(v$bumps) # Bumps sd
SNR <- 2 # Fix our SNR
sigma <- ssig/SNR
eps <- arima.sim(n=1024, model=list(ar=0.99, ma=1))
eps <- sigma*eps/sqrt(var(eps))
y <- v$bumps + eps
plot(x, y, type = "l", ylab = "Correlated Noise bumps")

ywd <- wd(y)# Take wavelet transform #
#
# Threshold BY LEVEL, but return the threshold
# value so we can print it out #
ywdT <- threshold(ywd, by.level=TRUE, policy="universal",return.thresh=TRUE)
#
# Now actually apply the threshold and invert
#
ywr <- wr(threshold(ywd, by.level=TRUE, policy="universal")) #
# Plot the denoised version and the original
#
yl <- range(c(ywr, v$bumps))
plot(x, ywr, ylim=yl, type="l")
lines(x, v$bumps, lty=2)

#extra example (from WaveThresh that I saw before) EbayesThresh
eps <- arima.sim(n=1024, model=list(ar=0.25, ma=1))
n<-length(eps)
y<-doppler(1:n/n)+eps/10
ywd <- wd(y)
ywdT <- threshold(ywd, by.level=TRUE, policy="universal",return.thresh=TRUE)
ywr <- wr(threshold(ywd, by.level=TRUE, policy="universal")) #
yl <- range(c(y,ywr))
plot(1:n, y, ylim=yl, type="l",col=8,lwd=1)
lines(1:n, ywr, lty=1,col=2,lwd=2)
lines(doppler(1:n/n),col=4,lwd=2)
