#' Clean up OM input (checks lengths of vectors for time-varying input) for use in \code{\link{GeMS}}
#'
#' @param input Vector or value to be "cleaned"
#' @param SimYear Number of years in the simulation
#'
#' @return Vector of the same length or greater than SimYear
#'
#' @export
#' @examples
#' \dontrun{
#' CheckReferencePoints("Cod_Base_CTL",ReadCTLfile("Cod_Base_CTL.csv"),"~/GeMS/Examples/Cod_1_Production")
#' }

CleanInput<-function(input,SimYear)
{
if(length(input)==1)
 output		<-rep(input,SimYear)
if(length(input)>1) {
 output		<-input[1:SimYear]
}
if(length(output)<SimYear) {
 stop("Check lengths of input vectors match (or are greater than) number of years.")
}
return(output)
}