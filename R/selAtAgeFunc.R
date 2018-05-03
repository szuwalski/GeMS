#' Selectivity-at-age
#' 
#' Calculates selectivity given lengths and growth parameters for use in \code{\link{PullTimeVary}}
#' 
#' @param inputLen Vector of lengths
#' @param K Growth rate (von Bertalanffy growth equation)
#' @param Linf Asymptotic maximum length (von Bertalanffy growth equation)
#' @param t0 Theoretical age at length=0 (von Bertalanffy growth equation)
#' 
#' @return Values of selectivity at age
#'
#' @export
#' @examples
#' \dontrun{
#' run_GeMS(c("Cod_LowProd_CTL","Cod_Base_CTL","Cod_HighProd_CTL"), MSEdir="~/GeneralMSE/Examples/Cod_1_Production")
#' }
selAtAgeFunc<-function(inputLen,K,Linf,t0)
 {
  selAtAge<-((log(1-inputLen/Linf)/-K))+t0
 return(selAtAge)
 }