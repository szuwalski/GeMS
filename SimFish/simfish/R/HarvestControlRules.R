#' \code{HarvestControlRule} runs a harvest control rule
#'
#' @param FMSY Fmsy
#' @param BMSY Bmsy
#' @param ExploitBio the exploitable biomass
#' @param SpawnBio the spawning biomass
#' @param alpha alpha something
#' @param beta beta something
#' @param HarvestControl form of harvest control rule
#' @param ConstantCatch constant catch
#' @param ConstantF constant F
#'
#' @return a TAC
#' @export
HarvestControlRule<-function(FMSY,BMSY,ExploitBio,SpawnBio,alpha,beta,HarvestControl,ConstantCatch=0,ConstantF=0)
{
  if(HarvestControl==1)
    outTAC<-ExploitBio*(FMSY)
  if(HarvestControl==2)
    outTAC<-ConstantCatch
  if(HarvestControl==3)
    outTAC<-ExploitBio*(ConstantF)

  if(HarvestControl==4)
  {
    useF	<-FMSY
    if(SpawnBio<BMSY)
    {
      useF	<-0
      if(SpawnBio>(BMSY*beta))
        useF	<-FMSY-((FMSY/(BMSY-beta*BMSY))*(BMSY-SpawnBio))
    }
    outTAC<-ExploitBio*useF
  }
  return(outTAC)
}
