library("wavethresh")
library("EbayesThresh")

nason<-function(y){
    ywd <- wd(y)
    ywdT <- threshold(ywd, by.level=TRUE, policy="universal",return.thresh=TRUE)
    ywr <- wr(threshold(ywd, by.level=TRUE, policy="universal"))
    list(yhat=ywr)
}




