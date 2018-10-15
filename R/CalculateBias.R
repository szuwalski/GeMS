#' Calculate Bias
#' @description Calculates bias in spawning stock biomass estimates for use in \code{\link{AgeStructureComp}}
#'
#' @param RetroOuts Output from \code{\link{CheckRetro}}
#'
#' @return Matrix of biases
#'
#' @export
#' @examples
#' \dontrun{
#' 
#' }

CalculateBias<-function(RetroOuts)
{
 preds<-RetroOuts[[1]]
 BiasSSB<-matrix(ncol=dim(preds)[3],nrow=dim(preds)[1])
 for(x in 1:dim(preds)[3])
  {
   tempPreds<-preds[,,x]
   tempTrue<-RetroOuts[[2]][x,]
   if(dim(preds)[1]>1) 
  	for(y in 1:dim(preds)[1])
   		BiasSSB[y,x]<-(tempPreds[y,(ncol(tempPreds)-y)]-tempTrue[(ncol(tempPreds)-y)])/tempTrue[(ncol(tempPreds)-y)]
   if(dim(preds)[1]==1)
   	BiasSSB[,x]<-(tempPreds[(length(tempPreds)-1)]-tempTrue[(length(tempPreds)-1)])/tempTrue[(length(tempPreds)-1)]
   }
 return(BiasSSB)
}
