#' Calculates Mohn's Rho for use in \code{\link{AgeStructureComp}}
#'
#' @param RetroOuts Output from \code{\link{CheckRetro}}
#'
#' @return Matrix of Mohn's Rhos
#'
#' @export
#' @examples
#' \dontrun{
#' 
#' }

CalculateMohn<-function(RetroOuts)
{
 preds<-RetroOuts[[1]]
 MohnsRho<-matrix(ncol=dim(preds)[3],nrow=(dim(preds)[1]-1))
 for(x in 1:dim(preds)[3])
  {
   tempPreds<-preds[,,x]
   tempTrue<-preds[1,,x]
   if(dim(preds)[1]>1) {
  	for(y in 2:dim(preds)[1]) {
   		MohnsRho[y-1,x]<-(tempPreds[y,(ncol(tempPreds)-y)]-tempTrue[(ncol(tempPreds)-y)])/tempTrue[(ncol(tempPreds)-y)]
  	}
   }
   if(dim(preds)[1]==1) {
   	MohnsRho<-NA
   }
  }
 return(MohnsRho)
}
