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
yyirregwd <- irregwd(yygrid)
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
## Not run: plot(yygrid$gridt, yy.wr, type="l")
#
# And superimpose the original data!
#
## Not run: points(tt, yy)
#
# This is sort of \code{Doppler} like!

#WaveThresh

#######################################
n=10^2
tt <- sort(runif(n))
#tt <- seq(0.001953125,0.9980469,length.out=n)
mu<-doppler(tt)
yy <-mu + arima.sim(n=100,list(ar=c(0.5)),sd=0.35)
#yy <-mu + rnorm(n, 0, 0.35)

yygrid <- makegrid(t=tt, y=yy)
yyirregwd <- irregwd(yygrid)
yy.thresh.sure <- threshold(yyirregwd, policy="sure", type="soft", dev=madmad) # policy sure is suppose to deal with correlation

yy.wr <- wr(yy.thresh.sure)
plot(yygrid$gridt, yy.wr, type="l")
points(tt, yy)
lines(tt,mu,col=4,lwd=3)

#plot(arima.sim(n=100,list(ar=c(0.5)),sd=0.35));lines(rnorm(n, 0, 0.35),col=2)
