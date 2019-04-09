#' Production Model
#' Fits a production model to input data
#' 
#' @param x Vector of production model parameters, c(carrying capacity, population growth rate, initial biomass as a proportion of carrying capacity)
#' @param CatchData Vector of catch
#' @param IndexData Vector of abundance indices
#' @param EstInit 1/0 Estimate initial biomass as proportion 
#' @param InitShape Initial value of shape parameter for Pella-Tomlinson, 1=Schaefer; 0=Fox
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
ProdMod<-function(x,CatchData,IndexData,estInit=0,InitShape=1)
{
	x<-exp(x)
	if(sum(is.na(x))>0 | sum(x<0)>0) return(Inf)

	predBio<-ProdModPlot(x,CatchData,IndexData,plots=0,estInit=estInit,InitShape=InitShape)[1:length(IndexData)]
	
	SSQ<-sum((predBio-IndexData)^2)
	return(SSQ)
}

