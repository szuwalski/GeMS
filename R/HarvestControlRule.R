#' Harvest Control Rule
#' Used within the \code{\link{GeMS}} function.
#' 
#' @param FMSY Fishing mortality at MSY
#' @param BMSY Biomass at MSY
#' @param ExploitBio Exploitable biomass
#' @param SpawnBio Spawning biomass
#' @param alpha Currently not implemented
#' @param beta Slope of harvest control rule
#' @param HarvestControl Harvest control rule (options in CTL file)
#' @param ConstantCatch Currently not implemented
#' @param ConstantF Currently not implemented

#' @param RetroPeels Number of years to peel for retrospective analysis
#' @param CTLNameList Vector of CTL file names
#' @param MSEdir Directory containing CTL files
#' @param plotNames Vector of the same length as CTLNameList with "prettier" names for plotting
#'
#' @return Value of Total Allowable Catch
#'
#' @export
#' 
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
