#' Harvest Control Rule
#' Used within the \code{\link{GeMS}} function.
#' 
#' @param FMSY Fishing mortality at MSY
#' @param BMSY Biomass at MSY
#' @param ExploitBio Exploitable biomass
#' @param SpawnBio Spawning biomass
#' @param alpha Slope of harvest control rule
#' @param beta Fraction of BMSY below which fishing mortality is set to zero
#' @param HarvestControl Harvest control rule (User-input HCR; 1=FMSY; 2=Constant Catch; 3=Constant F; 4=Sloped HCR)
#' @param ConstantCatch User-input catch (option in CTL file)
#' @param ConstantF User-input F (option in CTL file)

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
   if(SpawnBio<(BMSY*beta))
    {
    useF	<-0
    if(SpawnBio>(BMSY*beta)) 
      useF	<-FMSY-((FMSY/(BMSY-alpha*BMSY))*(BMSY-SpawnBio))
     }
   outTAC<-ExploitBio*useF
  }
 return(outTAC)
}
