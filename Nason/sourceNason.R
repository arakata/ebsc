library("wavethresh")
library("EbayesThresh")

#nason<-function(y){
#    ywd <- wd(y)
#    ywdT <- threshold(ywd, by.level=TRUE, policy="universal",return.thresh=TRUE)
#    ywr <- wr(threshold(ywd, by.level=TRUE, policy="universal"))
#    list(yhat=ywr)
#}

nason<-function(y){    
    t<-seq(0,1,length.out=length(y))
    data<-makegrid(t=t,y=y)
    y<-data$gridy
    ywd <- wd(y)
    ywdT <- threshold(ywd, by.level=TRUE, policy="universal",return.thresh=TRUE)
    ywr <- wr(threshold(ywd, by.level=TRUE, policy="universal"))
    list(t=t,yhat=ywr)
}





