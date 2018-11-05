#' Production Model
#' Fits a production model to input data
#' 
#' @param x Vector of production model parameters, c(carrying capacity, population growth rate, initial biomass as a proportion of carrying capacity)
#' @param CatchData Vector of catch
#' @param IndexData Vector of abundance indices
#' @param EstInit 1/0 Estimate initial biomass as proportion 
#' 
#' @return Sum Squared Errors from fitting production model
#'
#' @export
#' @examples
#' \dontrun{
#' K <- 10000
#' r <- 0.4
#' # Not needed if not estimating initial biomass
#' InitBio <- 0.95
#' x <- c(K,r,InitBio)
#' Catch <- c(500, 900, 750)
#' Index <- c(0.5, 0.9, 0.75)
#' ProdMod(x, Catch, Index, estInit=1)
#' }
#' 
ProdMod<-function(x,CatchData,IndexData,estInit=0)
{
	x<-exp(x)
	K<-x[1]
	r<-x[2]
	if(sum(is.na(x))>0 | sum(x<0)>0) return(Inf)
		predBio<-rep(0,length(IndexData))
		predBio[1]<-K
		if(estInit==1)
		  predBio[1]<-K*x[3]
		for(i in 2:length(CatchData))
		{
 			if(predBio[i-1]<=0) predBio[i] <- 0
 			if(predBio[i-1]>0) {
 			 predBio[i]<-predBio[i-1]+r*predBio[i-1]*(1-predBio[i-1]/K)-CatchData[i]
 			 if(predBio[i]<0) predBio[i] <- 0
 			}
		}

		SSQ<-sum((predBio-IndexData)^2)
		return(SSQ)
}