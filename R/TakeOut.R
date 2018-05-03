#' Extract values from given string for use in \code{\link{ReadCTLfile}}. Uses \code{grep} function.
#'  
#' @param SearchTerm Term to be searched
#' @param SearchPool String to be searched over
#' 
#' @return Vector of numbers or characters.
#'
#' @export
#' @examples
#' \dontrun{
#' TakeOut(SearchTerm="Nsim",SearchPool=as.matrix(read.csv(OMName)))
#' }

TakeOut<-function(SearchTerm,SearchPool)
{
 TakeIndex<-grep(SearchTerm,SearchPool[,3])[1]
 temp<-suppressWarnings(as.numeric(SearchPool[TakeIndex,1]))
 if(is.na(temp))
  temp<-eval(parse(text=SearchPool[TakeIndex,1]))
 return(temp)
}