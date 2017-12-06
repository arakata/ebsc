COiterations<-function(num.knots,degree,y){
	ndx=num.knots+1
	dx=(max(y)-min(y))/ndx
	knots=seq(min(y)-degree*dx,max(y)+degree*dx,by=dx)
	B=spline.des(knots,y,degree+1,0*y,outer.ok=T)$design
	y.hat=B%*%solve(t(B)%*%B)%*%t(B)%*%y
	return(y.hat=y.hat)
}