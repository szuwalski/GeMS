#' \code{CleanInput} extends parameters for the right amount of time
#'
#' @param input the thing to be extended
#' @param SimYear the number of years to extend
#'
#' @return the extended object
#' @export
CleanInput<-function(input,SimYear)
{
  if(length(input)==1)
    output		<-rep(input,SimYear)
  if(length(input)>1)
    output		<-input
  return(output)
}