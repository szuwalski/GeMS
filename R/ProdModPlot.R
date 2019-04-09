#' Production Model Plot
#' 
#' @param x Vector of production model parameters, c(carrying capacity, population growth rate, initial biomass as a proportion of carrying capacity)
#' @param CatchData Vector of catch
#' @param IndexData Vector of abundance indices
#' @param plots 1/0 To plot, or not to plot
#' @param EstInit 1/0 Estimate initial biomass as proportion 
#' 
#' @return Plot and vector of predicted biomass.
#'
#' @export
#' @examples
#' \dontrun{
#' K <- log(10000)
#' r <- 0.4
#' # Not needed if not estimating initial biomass
#' InitBio <- 0.95
#' x <- c(K,r,InitBio)
#' Catch <- c(500, 900, 750)
#' Index <- c(0.5, 0.9, 0.75)
#' ProdModPlot(x, Catch, Index, estInit=1)
#' }
#' 
ProdModPlot<-function(x,CatchData,IndexData,plots=0,estInit=0,InitShape=1)
{

	K<-x[1]
	r<-x[2]

	predBio<-rep(0,length(IndexData)+1)
	predBio[1]<-K

	if(!InitShape%in%c(1,0))
		p<-x[3]
	if(InitShape%in%c(1,0))
		p<-InitShape

	if(p==0) p<-1E-10

	if(estInit==1 & InitShape%in%c(1,0))
		predBio[1]<-K*x[3]
	if(estInit==1 & !InitShape%in%c(1,0))
		predBio[1]<-K*x[4]


	for(i in 2:length(CatchData))
	{
 		if(predBio[i-1]<=0) predBio[i] <- 0
 		if(predBio[i-1]>0) {
 		 predBio[i]<-predBio[i-1]+(r/p)*predBio[i-1]*(1-(predBio[i-1]/K)^p)-CatchData[i]
 		 if(predBio[i]<0) predBio[i] <- 0
 		}
	}

	if(plots==1)
	{
		plot(IndexData[1:(length(IndexData)-1)],ylim=c(0,max(IndexData,na.rm=T)),las=1,ylab="")
		lines(predBio[1:(length(predBio)-1)])
		lines(CatchData[1:(length(CatchData)-1)],lty=2,col=2)
	}
	return(predBio)
}