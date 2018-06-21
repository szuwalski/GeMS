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
#' MSEdir <- "~/GeneralMSE/Examples/Cod_5_AgeStructure"
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
TRU			<-readLines(file.path(MSEdir,IndSimFolder,"TrueQuantities.DAT"))

#==true
temp		<-grep("fishing mortality",TRU)
TrueFmort[,,p]<-matrix(suppressWarnings(as.numeric(unlist(strsplit(TRU[(temp+1):(temp+Nsim)],split=" ")))),ncol=SimYear,byrow=T)

#==true Catch
temp		<-grep("Catch",TRU)
TrueCatch[,,p]<-matrix(suppressWarnings(as.numeric(unlist(strsplit(TRU[(temp+1):(temp+Nsim)],split=" ")))),ncol=SimYear,byrow=T)

#==true Spawning biomass
temp		<-grep("spawning biomass",TRU)
TrueSpbio[,,p]<-matrix(suppressWarnings(as.numeric(unlist(strsplit(TRU[(temp+1):(temp+Nsim)],split=" ")))),ncol=SimYear,byrow=T)

temp		<-grep("recruitment",TRU)
TrueRec[,,p]<-matrix(suppressWarnings(as.numeric(unlist(strsplit(TRU[(temp+1):(temp+Nsim)],split=" ")))),ncol=SimYear,byrow=T)

#==OFL 
temp			<-grep("OFL",TRU)
TrueOFL[,,p]	<-matrix(suppressWarnings(as.numeric(unlist(strsplit(TRU[(temp+1):(temp+Nsim)],split=" ")))),ncol=SimYear,byrow=T)

for(n in 1:Nsim)
{
IndSimFolder	<-file.path(DrawDir,n,SimYear)
PAR			<-readLines(file.path(MSEdir,IndSimFolder,"simass.par"))

#==fishing mortality
temp		<-grep("log_avg_fmort_dir",PAR)[1]
AvgF		<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))
temp		<-grep("fmort_dir_dev",PAR)[1]
fDevs		<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))[-1]
FishMort[n,AssYear:(SimYear-1),p]<-exp(AvgF+fDevs)

#==est spawning biomass
IndSimFolder2			<-file.path(DrawDir,n,SimYear)
pullREP				<-readLines(file.path(MSEdir,IndSimFolder2,"simass.rep"))
temp					<-grep("spawning biomass",pullREP)
EstSpbio[n,AssYear:(SimYear-1),p]			<-as.numeric((unlist(strsplit(pullREP[temp+1],split=" "))))[2:(SimYear-AssYear+1)]

#==recruitment
temp		<-grep("mean_log_rec",PAR)[1]
AvgRec	<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))
temp		<-grep("rec_dev",PAR)[1]
RecDevs	<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))[-1]
Recruitment[n,AssYear:(SimYear-1),p]<-exp(AvgRec+RecDevs)

for(w in (InitYear+1):SimYear)
{
IndSimFolder2			<-file.path(DrawDir,n,w)
pullREP				<-readLines(file.path(MSEdir,IndSimFolder2,"simass.rep"))
temp					<-grep("OFL",pullREP)[2]
EstOFL[n,w,p]			<-as.numeric((unlist(strsplit(pullREP[temp+1],split=" "))))
}

#==timevary
#==growth combos
#==write a function already that pulls the variable...

if(out$OM$TimeVaryLinf>0 & out$OM$TimeVaryGrowthK <=0)
{
temp			<-grep("log_avg_Linf",PAR)[1]
log_avg_Linf	<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))
temp			<-grep("Linf_dev",PAR)[1]
Linf_dev		<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))[-1]
Linf			<-exp(log_avg_Linf+Linf_dev)

GrowthK		<-rep(out$OM$VonKn,SimYear)
}

if(out$OM$TimeVaryLinf<=0 & out$OM$TimeVaryGrowthK >0)
{
temp			<-grep("log_avg_GrowthK",PAR)[1]
log_avg_GrowthK	<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))
temp			<-grep("GrowthK_dev",PAR)[1]
GrowthK_dev		<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))[-1]
GrowthK		<-exp(log_avg_GrowthK+GrowthK_dev)

Linf			<-rep(out$OM$LinfN,SimYear)
}

if(out$OM$TimeVaryLinf>0 & out$OM$TimeVaryGrowthK>0)
{
temp			<-grep("log_avg_Linf",PAR)[1]
log_avg_Linf	<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))
temp			<-grep("Linf_dev",PAR)[1]
Linf_dev		<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))[-1]
Linf			<-exp(log_avg_Linf+Linf_dev)

temp			<-grep("log_avg_GrowthK",PAR)[1]
log_avg_GrowthK	<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))
temp			<-grep("GrowthK_dev",PAR)[1]
GrowthK_dev		<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))[-1]
GrowthK		<-exp(log_avg_GrowthK+GrowthK_dev)
}

if(out$OM$TimeVaryLinf<=0 & out$OM$TimeVaryGrowthK<=0)
{
Linf			<-rep(out$OM$LinfN,SimYear)
GrowthK		<-rep(out$OM$VonKn,SimYear)
}

for(x in 1:SimYear)
 GrowthVary[x,,n,p]<-Linf[x]*(1-exp(-GrowthK[x]*(Ages-t0)))

#==natural mortality
if(out$OM$TimeVaryM >0)
{
temp			<-grep("log_avg_NatM",PAR)[1]
avgM			<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))
temp			<-grep("m_dev",PAR)[1]
mDevs			<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))[-1]
NatMvary[n,out$OM$start_assessment:(SimYear-1),p]	<-exp(avgM+mDevs)
}

if(out$OM$TimeVaryM <=0 )
{
  temp			<-grep("log_avg_NatM",PAR)[1]
  avgM			<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))
  NatMvary[n,out$OM$start_assessment:(SimYear-1),p]	<-exp(avgM)
}

#==fishing selectivity combos
# selAtAgeFunc(sel50n[1],VonKn[1],LinfN[1],t0n[1])
# selAtAgeFunc(sel95n[1],VonKn[1],LinfN[1],t0n[1])


if(out$OM$TimeVarySel50 >0 & out$OM$TimeVarySel95 <=0)
{
temp			<-grep("SelPars50",PAR)[1]
SelPars50		<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))
temp			<-grep("SelPars_dev50",PAR)[1]
SelPars_dev50	<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))[-1]
Sel50			<-SelPars50+SelPars_dev50

temp			<-grep("SelPars95",PAR)[1]
SelPars95		<-rep(as.numeric(unlist(strsplit(PAR[temp+1],split=" "))),SimYear)
}

if(out$OM$TimeVarySel50 <=0 & out$OM$TimeVarySel95 >0)
{
temp			<-grep("SelPars95",PAR)[1]
SelPars95		<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))
temp			<-grep("SelPars_dev95",PAR)[1]
SelPars_dev95	<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))[-1]
Sel95			<-SelPars95+SelPars_dev95

temp			<-grep("SelPars50",PAR)[1]
SelPars50		<-rep(as.numeric(unlist(strsplit(PAR[temp+1],split=" "))),SimYear)
}

if(out$OM$TimeVarySel50>0 & out$OM$TimeVarySel95 >0)
{
temp			<-grep("SelPars95",PAR)[1]
SelPars95		<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))
temp			<-grep("SelPars_dev95",PAR)[1]
SelPars_dev95	<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))[-1]
Sel95			<-SelPars95+SelPars_dev95

temp			<-grep("SelPars50",PAR)[1]
SelPars50		<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))
temp			<-grep("SelPars_dev50",PAR)[1]
SelPars_dev50	<-as.numeric((unlist(strsplit(PAR[temp+1],split=" "))))[-1]
Sel50			<-SelPars50+SelPars_dev50
}

if(out$OM$TimeVarySel50<=0 & out$OM$TimeVarySel95 <=0)
{
temp			<-grep("SelPars50",PAR)[1]
Sel50		<-rep(as.numeric(unlist(strsplit(PAR[temp+1],split=" "))),SimYear)
temp			<-grep("SelPars95",PAR)[1]
Sel95		<-rep(as.numeric(unlist(strsplit(PAR[temp+1],split=" "))),SimYear)
}

for(x in 1:SimYear)
 SelVary[x,,n,p]	<-1/( 1 + exp( -1*log(19)*(Ages-Sel50[x])/(Sel95[x]-Sel50[x])))

}
#==============================================
# Pull the 'true' values for each process
#==============================================
tGrowthK		<-out$OM$VonKn
if(length(tGrowthK)==1)
 tGrowthK		<-rep(tGrowthK,SimYear)
tLinf			<-out$OM$LinfN
if(length(tLinf)==1)
 tLinf		<-rep(tLinf,SimYear)

for(x in 1:SimYear)
 TrueGrow[x,,p]	<-tLinf[x]*(1-exp(-tGrowthK[x]*(Ages-t0)))
	
 sel50n<-CleanInput(out$OM$sel50n,out$OM$SimYear)
 sel95n<-CleanInput(out$OM$sel95n,out$OM$SimYear)
 vonKn<-CleanInput(out$OM$vonKn,out$OM$SimYear)
 LinfN<-CleanInput(out$OM$LinfN,out$OM$SimYear)
 t0n<-CleanInput(out$OM$t0n,out$OM$SimYear)

 tSel50<-selAtAgeFunc(sel50n,VonKn,LinfN,t0n)
 tSel95<-selAtAgeFunc(sel95n,VonKn,LinfN,t0n)

#if(length(tSel50)==1)
# tSel50		<-rep(tSel50,out$OM$SimYear)
#if(length(tSel95)==1)
# tSel95		<-rep(tSel95,out$OM$SimYear)

for(x in 1:SimYear)
 TrueSel[x,,p]	<-1/( 1 + exp( -1*log(19)*(Ages-tSel50[x])/(tSel95[x]-tSel50[x])))

TrueM			<-CleanInput(out$OM$NatMn,out$OM$SimYear)
#if(length(TrueM)==1)
# TrueM		<-rep(TrueM,SimYear)



#==pull observed and fitted values
for(n in 1:Nsim)
{
  IndSimFolder	<-file.path(DrawDir,n,SimYear)
  REP			<-readLines(file.path(MSEdir,IndSimFolder,"simass.rep"))
  DAT			<-readLines(file.path(MSEdir,IndSimFolder,"simass.dat"))
  
temp<-grep("survey years",DAT)[1]
yearsDat[p]<-as.numeric(unlist(strsplit(DAT[temp+1],split=" ")))

temp<-grep("pred survey bio",REP)
predSurvBio[n,,p]<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat[p]+1)]
temp<-grep("obs survey bio",REP)
obsSurvBio[n,,p]<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat[p]+1)]

temp<-grep("pred cpue bio",REP)
predCpueBio[n,,p] <-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat[p]+1)]
temp<-grep("obs cpue bio",REP)
cpueIndex[n,,p] <-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat[p]+1)]

temp<-grep("pred catch bio",REP)
predCatchBio[n,,p]<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat[p]+1)]
temp<-grep("obs catch bio",REP)
obsCatchBio[n,,p]<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat[p]+1)]

}

}

list(Recruitment,TrueRec,FishMort,TrueFmort,GrowthVary,TrueGrow,NatMvary,TrueM,SelVary,
     TrueSel,EstOFL,TrueOFL,TrueCatch,TrueSpbio,EstSpbio,
     predSurvBio,obsSurvBio,predCpueBio,cpueIndex,predCatchBio,obsCatchBio)
}
