#' Used to pull parameters from age-structured estimation model output. For use in \code{\link{AgeStructureComp}}
#' 
#' @param out List output from \code{\link{ReadCTLfile}}
#' @param MSEdir Directory containing CTL files
#' @param CTLNameList Vector of CTL names
#'
#' @return List of quantities
#'
#' @export
#' @examples
#' \dontrun{
#' MSEdir <- "~/GeMS/inst/extdata/Cod_5_AgeStructure"
#' CTLName <- "Cod_Age_Structure_CTL"
#' out <- ReadCTLfile(CTLName)
#' PullTimevary(out=out,
#'              MSEdir=MSEdir,
#'              CTLNameList=CTLName)
#' }
#' 
PullTimevary<-function(out,MSEdir,CTLNameList)
{
Nsim			<-out$OM$Nsim			# number of simulations to do in the MSE
SimYear		<-out$OM$SimYear			# total number of years in simulation
InitYear		<-out$OM$InitYear			# year in which MSE starts (i.e. the number of years of data available)
MaxAge		<-out$OM$MaxAge
Ages			<-seq(1,MaxAge)
t0			<-out$OM$t0s		# check this later when going to spatial
AssYear <- out$OM$start_assessment

Recruitment	<-array(dim=c(nrow=(Nsim),ncol=SimYear-1,length(CTLNameList)))
FishMort	<-array(dim=c(nrow=(Nsim),ncol=SimYear-1,length(CTLNameList)))
TrueRec	<-array(dim=c(nrow=(Nsim),ncol=SimYear,length(CTLNameList)))
TrueFmort	<-array(dim=c(nrow=(Nsim),ncol=SimYear,length(CTLNameList)))
TrueCatch	<-array(dim=c(nrow=(Nsim),ncol=SimYear,length(CTLNameList)))
TrueSpbio	<-array(dim=c(nrow=(Nsim),ncol=SimYear,length(CTLNameList)))
EstSpbio	<-array(dim=c(nrow=(Nsim),ncol=SimYear,length(CTLNameList)))
NatMvary	<-array(dim=c(nrow=(Nsim),ncol=SimYear-1,length(CTLNameList)))
SelVary	<-array(dim=c(SimYear,MaxAge,Nsim,length(CTLNameList))) 
GrowthVary	<-array(dim=c(SimYear,MaxAge,Nsim,length(CTLNameList)))
TrueGrow	<-array(dim=c(nrow=SimYear,ncol=MaxAge,length(CTLNameList)))
TrueSel	<-array(dim=c(nrow=SimYear,ncol=MaxAge,length(CTLNameList)))
TrueM   <-array(dim=c(nrow=(Nsim),ncol=SimYear,length(CTLNameList)))
TrueOFL	<-array(dim=c(nrow=(Nsim),ncol=SimYear,length(CTLNameList)))
EstOFL	<-array(dim=c(nrow=(Nsim),ncol=SimYear,length(CTLNameList)))

predSurvBio<-array(dim=c(nrow=(Nsim),ncol=length(AssYear:SimYear)-1,length(CTLNameList)))
obsSurvBio<-array(dim=c(nrow=(Nsim),ncol=length(AssYear:SimYear)-1,length(CTLNameList)))
predCpueBio<-array(dim=c(nrow=(Nsim),ncol=length(AssYear:SimYear)-1,length(CTLNameList)))
cpueIndex<-array(dim=c(nrow=(Nsim),ncol=length(AssYear:SimYear)-1,length(CTLNameList)))
predCatchBio<-array(dim=c(nrow=(Nsim),ncol=length(AssYear:SimYear)-1,length(CTLNameList)))
obsCatchBio<-array(dim=c(nrow=(Nsim),ncol=length(AssYear:SimYear)-1,length(CTLNameList)))
yearsDat<-rep(0,length(CTLNameList))


for(p in seq_along(CTLNameList))
{
DrawDir		<-CTLNameList[p]
out			<-ReadCTLfile(file.path(MSEdir,paste0(CTLNameList[p])))

IndSimFolder	<-file.path(DrawDir,Nsim,SimYear)
TRU			<-ReadTruDat(IndSimFolder)

#==true Fmort
TrueFmort[,,p]<-TRU$fishing_mortality
  
#==true Catch
TrueCatch[,,p]<-TRU$Catch

#==true Spawning biomass
TrueSpbio[,,p]<-TRU$spawning_biomass
  
TrueRec[,,p]<-TRU$recruitment
  
#==OFL 
TrueOFL[,,p]	<-TRU$OFL
  
for(n in 1:Nsim)
{
IndSimFolder	<-file.path(DrawDir,n,SimYear)
PAR			<-ReadPar(IndSimFolder)

#==fishing mortality
AvgF		<- PAR$log_avg_fmort_dir
fDevs		<- PAR$fmort_dir_dev
FishMort[n,AssYear:(SimYear-1),p]<-exp(AvgF+fDevs)

#==est spawning biomass
IndSimFolder2			<-file.path(DrawDir,n,SimYear)
REP				    <-ReadRep(file.path(MSEdir,IndSimFolder2))
EstSpbio[n,AssYear:(SimYear-1),p]			<-REP$spawning_biomass

#==recruitment
AvgRec	<-PAR$mean_log_rec
RecDevs	<-PAR$rec_dev
Recruitment[n,AssYear:(SimYear-1),p]<-exp(AvgRec+RecDevs)

for(w in (InitYear+1):SimYear)
{
IndSimFolder2			<-file.path(DrawDir,n,w)
REP				        <-ReadRep(file.path(MSEdir,IndSimFolder2))
EstOFL[n,w,p]			<-REP$OFL
}

#==timevary
#==growth combos
#==write a function already that pulls the variable...
log_avg_Linf	<-PAR$log_avg_Linf
Linf_dev		<-PAR$Linf_dev
Linf			<-exp(log_avg_Linf+Linf_dev)

log_avg_GrowthK	<-PAR$log_avg_GrowthK
GrowthK_dev		<-PAR$GrowthK_dev
GrowthK		<-exp(log_avg_GrowthK+GrowthK_dev)

for(x in 1:SimYear)
 GrowthVary[x,,n,p]<-Linf[x]*(1-exp(-GrowthK[x]*(Ages-t0)))

#==natural mortality
avgM			<-PAR$log_avg_NatM
mDevs			<-PAR$m_dev
NatMvary[n,out$OM$start_assessment:(SimYear-1),p]	<-exp(avgM+mDevs)

#==fishing selectivity combos
# selAtAgeFunc(sel50n[1],VonKn[1],LinfN[1],t0n[1])
# selAtAgeFunc(sel95n[1],VonKn[1],LinfN[1],t0n[1])

SelPars50		<- PAR$SelPars50
SelPars_dev50	<- PAR$SelPars_dev50
Sel50			<-SelPars50+SelPars_dev50

SelPars95		<-PAR$SelPars95
SelPars_dev95	<-PAR$SelPars_dev95
Sel95			<-SelPars95+SelPars_dev95

for(x in 1:SimYear)
 SelVary[x,,n,p]	<-1/( 1 + exp( -1*log(19)*(Ages-Sel50[x])/(Sel95[x]-Sel50[x])))

}

#==============================================
# Pull the 'true' values for each process
#==============================================
tGrowthK		<-CleanInput(out$OM$VonKn, SimYear)
# if(length(tGrowthK)==1)
#  tGrowthK		<-rep(tGrowthK,SimYear)
tLinf			<-CleanInput(out$OM$LinfN,SimYear)
# if(length(tLinf)==1)
#  tLinf		<-rep(tLinf,SimYear)

for(x in 1:SimYear)
 TrueGrow[x,,p]	<-tLinf[x]*(1-exp(-tGrowthK[x]*(Ages-t0)))
	
 sel50n<-CleanInput(out$OM$sel50n,SimYear)
 sel95n<-CleanInput(out$OM$sel95n,SimYear)
 t0n<-CleanInput(out$OM$t0n,SimYear)

 tSel50<-selAtAgeFunc(sel50n,tGrowthK,tLinf,t0n)
 tSel95<-selAtAgeFunc(sel95n,tGrowthK,tLinf,t0n)

#if(length(tSel50)==1)
# tSel50		<-rep(tSel50,out$OM$SimYear)
#if(length(tSel95)==1)
# tSel95		<-rep(tSel95,out$OM$SimYear)

for(x in 1:SimYear)
 TrueSel[x,,p]	<-1/( 1 + exp( -1*log(19)*(Ages-tSel50[x])/(tSel95[x]-tSel50[x])))

TrueM			<-CleanInput(out$OM$NatMn,SimYear)
#if(length(TrueM)==1)
# TrueM		<-rep(TrueM,SimYear)



#==pull observed and fitted values
for(n in 1:Nsim)
{
  IndSimFolder	<-file.path(DrawDir,n,SimYear)
  REP			<-ReadRep(IndSimFolder)
  
  yearsDat[p]       <- REP$yearsDat
  predSurvBio[n,,p] <- REP$pred_survey_bio
  obsSurvBio[n,,p]  <- REP$obs_survey_bio
  predCpueBio[n,,p] <- REP$pred_cpue_bio
  cpueIndex[n,,p]   <- REP$obs_cpue_bio
  predCatchBio[n,,p]<- REP$pred_catch_bio
  obsCatchBio[n,,p] <- REP$obs_catch_bio
}

}

list(Recruitment=Recruitment,TrueRec=TrueRec,FishMort=FishMort,TrueFmort=TrueFmort,GrowthVary=GrowthVary,TrueGrow=TrueGrow,NatMvary=NatMvary,TrueM=TrueM,SelVary=SelVary,
     TrueSel=TrueSel,EstOFL=EstOFL,TrueOFL=TrueOFL,TrueCatch=TrueCatch,TrueSpbio=TrueSpbio,EstSpbio=EstSpbio,
     predSurvBio=predSurvBio,obsSurvBio=obsSurvBio,predCpueBio=predCpueBio,cpueIndex=cpueIndex,predCatchBio=predCatchBio,obsCatchBio=obsCatchBio)
}
