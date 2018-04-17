###################################################################
# Generalized spatial management strategy evaluation
# Two box, populations processes in each box can
# vary over time, movement between boxes can vary over 
# time, effort in boxes can be allocated asymmetrically
# production and age-structured estimation models can be tested
# written by Cody Szuwalski, 11/2015
###################################################################

GeMS<-function(out,CreateFolderName,MSEdir,GeMSdir=NULL,silent=F,ADsilent=T,ADoptions=NA)
{

echofile <- file.path(MSEdir,paste0(CreateFolderName,"_echo.txt"))
if(!file.exists(echofile)) file.create(echofile)
cat(paste0("# Run initiated: ",Sys.time(),"\n"),file=echofile,append=F)
cat(paste0("# Best viewed in Excel or some tab-delimiting editor.","\n"),file=echofile,append=T)
cat("## Inputs and Derived Quantities \n",file=echofile,append=T)
cat(CreateFolderName,"\t # CreateFolderName \n",file=echofile,append=T)

#===============
#==Checking OS==
#===============
os <- .Platform$OS.type

if(os != "windows") {
  SimAssExec <- "SimAss"
  SimAssComm <- "./SimAss -nohess"
}

if(os == "windows") {
  SimAssExec <- "simass.exe"
  SimAssComm <- "simass -nohess"
}

if(ADsilent) {
  SimAssComm <- paste(SimAssComm, ">console.log", sep = " ") # saves console output to a file "console.log"
}

if(is.na(ADoptions!=NA)) {
  SimAssComm <- paste(SimAssComm, ADoptions, sep = " ") # ADMB options such as "-gbs" or "-cbs" (Memory Management)
}

#=======================
#==simulation controls==
#=======================
Nsim			<-out$OM$Nsim			# number of simulations to do in the MSE
cat(Nsim,"\t # Nsim \n",file=echofile,append=T)
SimYear		<-out$OM$SimYear			# total number of years in simulation
cat(SimYear,"\t # SimYear \n",file=echofile,append=T)
InitYear		<-out$OM$InitYear			# year in which MSE starts (i.e. the number of years of data available for initial assessment)
cat(InitYear,"\t # InitYear \n",file=echofile,append=T)
AssessmentType	<-out$OM$AssessmentType		# 0 = projection only, 1 = production, 2= agestructured
cat(AssessmentType,"\t # AssessmentType \n",file=echofile,append=T)
FisheryIndepenDat <-out$OM$FisheryIndepenDat
cat(FisheryIndepenDat,"\t # FisheryIndepenDat \n",file=echofile,append=T)
LifeHistoryPlots	<-out$OM$LifeHistoryPlots
EstimationPlots	<-out$OM$EstimationPlots	# this plots the diagnostics for each run of the estimation model (will be a lot of plots!)
AssessmentData	<-out$OM$AssessmentData		# this plots the diagnostics for each run of the estimation model (will be a lot of plots!)
PlotYieldCurve	<-out$OM$PlotYieldCurve
TwoPop		<-out$OM$TwoPop

PlotFolder<-file.path(MSEdir, "plots")
if(!file.exists(PlotFolder)) dir.create(PlotFolder)

#==sampling uncertainty 
CatchCVn	<-CleanInput(out$OM$CatchCVn,SimYear)
cat(CatchCVn,"\t # CatchCVn \n",file=echofile,append=T)
CatchCVs	<-CatchCVn
if(TwoPop>0) 
 CatchCVs	<-CleanInput(out$OM$CatchCVs,SimYear)

IndexCVn	<-CleanInput(out$OM$IndexCVn,SimYear)
cat(IndexCVn,"\t # IndexCVn \n",file=echofile,append=T)
IndexCVs	<-IndexCVn
if(TwoPop>0) 
 IndexCVs	<-CleanInput(out$OM$IndexCVs,SimYear)

LenSampleN	<-CleanInput(out$OM$LenSampleN,SimYear)
cat(LenSampleN,"\t # LenSampleN \n",file=echofile,append=T)
LenSampleS	<-LenSampleN
if(TwoPop>0) 
 LenSampleS	<-CleanInput(out$OM$LenSampleS,SimYear)

GrowthSDn	<-CleanInput(out$OM$GrowthSDn,SimYear)
cat(GrowthSDn,"\t # GrowthSDn \n",file=echofile,append=T)
GrowthSDs	<-GrowthSDn
if(TwoPop>0) 
 GrowthSDs	<-CleanInput(out$OM$GrowthSDs,SimYear)

#==========================================================
#=================population dynamics processes============
#==============================================+===========
MaxAge<-out$OM$MaxAge
cat(MaxAge,"\t # MaxAge \n",file=echofile,append=T)
Ages<-seq(1,MaxAge)

#==natural mortality================
 NatMn	<-CleanInput(out$OM$NatMn,SimYear)
 cat(NatMn,"\t # NatMn \n",file=echofile,append=T)
 NatMs	<-NatMn
 if(TwoPop>0) 
  NatMs	<-CleanInput(out$OM$NatMn,SimYear)

#==Length at age====================
VonKn		<-CleanInput(out$OM$VonKn,SimYear)
cat(VonKn,"\t # VonKn \n",file=echofile,append=T)
LinfN		<-CleanInput(out$OM$LinfN,SimYear)
cat(LinfN,"\t # LinfN \n",file=echofile,append=T)
t0n		<-CleanInput(out$OM$t0n,SimYear)
cat(t0n,"\t # t0n \n",file=echofile,append=T)
LenAtAgeN	<-matrix(nrow=SimYear,ncol=MaxAge)

LenAtAgeN[1,]<-LinfN[1]*(1-exp(-VonKn[1]*(Ages-t0n[1])))
LenAtAgeN[,1]<-LenAtAgeN[1,1]
for(i in 2:SimYear)
 for(j in 2:MaxAge)
  LenAtAgeN[i,j]<-LenAtAgeN[i-1,j-1]+(LinfN[i]-LenAtAgeN[i-1,j-1])*(1-exp(-VonKn[i]))
cat("# LenAtAgeN \n",file=echofile,append=T)
cat("# Ages", Ages,"\n",file=echofile,append=T)
write.table(row.names=F,col.names=F,LenAtAgeN,file=echofile,append=T)

LenAtAgeS<-LenAtAgeN

if(TwoPop>0)
{ 
VonKs		<-CleanInput(out$OM$VonKs,SimYear)
LinfS		<-CleanInput(out$OM$LinfS,SimYear)
t0s		<-CleanInput(out$OM$t0s,SimYear)
LenAtAgeS	<-matrix(nrow=SimYear,ncol=MaxAge)

LenAtAgeS[1,]<-LinfS[1]*(1-exp(-VonKs[1]*(Ages-t0s[1])))
LenAtAgeS[,1]<-LenAtAgeS[1,1]
for(i in 2:SimYear)
 for(j in 2:MaxAge)
  LenAtAgeS[i,j]<-LenAtAgeS[i-1,j-1]+(LinfS[i]-LenAtAgeS[i-1,j-1])*(1-exp(-VonKs[i]))
}

#==specify the number of length bins
 BinWidth		<-(max(LenAtAgeS,LenAtAgeN)*1.05)/(out$OM$LengthBinN)
 cat("\n",BinWidth,"\t # BinWidth \n",file=echofile,append=T)
 LengthBins		<-seq(0,max(LenAtAgeS,LenAtAgeN)*1.05,BinWidth)
 cat(LengthBins,"\t # LengthBins \n",file=echofile,append=T)
 LengthBinsMid	<-LengthBins[1:(length(LengthBins)-1)] + mean(LengthBins[1:2])
 LengthBinN		<-length(LengthBinsMid)
 cat(LengthBinN,"\t # LengthBinN \n",file=echofile,append=T)

#==maturity at age==========================
mat50n	<-CleanInput(out$OM$mat50n,SimYear)
cat(mat50n,"\t # mat50n \n",file=echofile,append=T)
mat95n	<-CleanInput(out$OM$mat95n,SimYear)
cat(mat95n,"\t # mat95n \n",file=echofile,append=T)
matureN	<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 matureN[i,]<-1/(1+exp(-1*log(19)*(Ages-mat50n[i])/(mat95n[i]-mat50n[i])))
cat("# matureN \n",file=echofile,append=T)
write.table(row.names=F,col.names=F,matureN,file=echofile,append=T)

matureS	<-matureN

if(TwoPop>0) 
{
mat50s	<-CleanInput(out$OM$mat50s,SimYear)
mat95s	<-CleanInput(out$OM$mat95s,SimYear)
matureS	<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 matureS[i,]<-1/(1+exp(-1*log(19)*(Ages-mat50s[i])/(mat95s[i]-mat50s[i])))
}

#=====================================
#==fishery selectivitiy============== 
sel50n	<-CleanInput(out$OM$sel50n,SimYear)
cat(sel50n,"\t # sel50n \n",file=echofile,append=T)
sel95n	<-CleanInput(out$OM$sel95n,SimYear)
cat(sel95n,"\t # sel95n \n",file=echofile,append=T)
vulnN		<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 vulnN[i,]	<-1/(1+exp(-1*log(19)*(LenAtAgeN[i,]-sel50n[i])/(sel95n[i]-sel50n[i])))

vulnN[vulnN<0.01]<-0
cat("\n # vulnN \n",file=echofile,append=T)
write.table(row.names=F,col.names=F,vulnN,file=echofile,append=T)
vulnS		<-vulnN
if(TwoPop>0) 
{
sel50s	<-CleanInput(out$OM$sel50s,SimYear)
sel95s	<-CleanInput(out$OM$sel95s,SimYear)
vulnS		<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 vulnS[i,]	<-1/(1+exp(-1*log(19)*(LenAtAgeN[i,]-sel50s[i])/(sel95s[i]-sel50s[i])))

vulnS[vulnS<0.01]<-0
}
#==index selectivity=====================
surv50n	<-CleanInput(out$OM$surv50n,SimYear)
cat(surv50n,"\t # surv50n \n",file=echofile,append=T)
surv95n	<-CleanInput(out$OM$surv95n,SimYear)
cat(surv95n,"\t # surv95n \n",file=echofile,append=T)
survSelN		<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 survSelN[i,]	<-1/(1+exp(-1*log(19)*(LenAtAgeN[i,]-surv50n[i])/(surv95n[i]-surv50n[i])))
cat("# survSelN \n",file=echofile,append=T)
write.table(row.names=F,col.names=F,survSelN,file=echofile,append=T)

survSelS	<-survSelN
surv50s	<-surv50n
surv95s	<-surv95n

if(TwoPop>0) 
{
surv50s	<-CleanInput(out$OM$surv50s,SimYear)
surv95s	<-CleanInput(out$OM$surv95s,SimYear)
survSelS	<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 survSelS[i,]	<-1/(1+exp(-1*log(19)*(LenAtAgeN[i,]-surv50s[i])/(surv95s[i]-surv50s[i])))
}

if(max(sel50n)>max(LinfN) | max(sel95n)>max(LinfN) | max(surv50n)>max(LinfN) | max(surv50n)>max(LinfN)) {stop("Gear is attempting to select fish larger than what's in the population (one or more of the selectivity parameters is greater than Linf)")}
if(TwoPop>0) if(max(sel50s)>max(LinfS) | max(sel95s)>max(LinfS) | max(surv50s)>max(LinfS) | max(surv50s)>max(LinfS)) {stop("Gear is attempting to select fish larger than what's in the population (one or more of the selectivity parameters is greater than Linf)")}


#==weight at age==========================
alphaN		<-CleanInput(out$OM$alphaN,SimYear)
cat(alphaN,"\t # alphaN \n",file=echofile,append=T)
betaN			<-CleanInput(out$OM$betaN,SimYear)
cat(betaN,"\t # betaN \n",file=echofile,append=T)
WeightAtAgeN	<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 WeightAtAgeN[i,]	<-alphaN[i]*LenAtAgeN[i,]^betaN[i]
cat("\n # WeightAtAgeN \n",file=echofile,append=T)
write.table(row.names=F,col.names=F,WeightAtAgeN,file=echofile,append=T)

WeightAtAgeS	<-WeightAtAgeN
if(TwoPop>0) 
{
alphaS		<-CleanInput(out$OM$alphaS,SimYear)
betaS			<-CleanInput(out$OM$betaS,SimYear)
WeightAtAgeS	<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 WeightAtAgeS[i,]	<-alphaS[i]*LenAtAgeS[i,]^betaS[i]
}

#==movement from box to box==
MaxMovingN	<-CleanInput(out$OM$MaxMovingN,SimYear)
cat(MaxMovingN,"\t # MaxMovingN \n",file=echofile,append=T)
Move50n	<-CleanInput(out$OM$Move50n,SimYear)
cat(Move50n,"\t # Move50n \n",file=echofile,append=T)
Move95n	<-CleanInput(out$OM$Move95n,SimYear)
cat(Move95n,"\t # Move95n \n",file=echofile,append=T)
MovementN	<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 MovementN[i,]	<-MaxMovingN[i]/(1+exp(-1*log(19)*(Ages-Move50n[i])/(Move95n[i]-Move50n[i])))
cat("\n # MovementN \n",file=echofile,append=T)
write.table(row.names=F,col.names=F,MovementN,file=echofile,append=T)

MovementS	<-MovementN
if(TwoPop>0) 
{
MaxMovingS	<-CleanInput(out$OM$MaxMovingS,SimYear)
Move50s	<-CleanInput(out$OM$Move50s,SimYear)
Move95s	<-CleanInput(out$OM$Move95s,SimYear)
MovementS	<-matrix(nrow=SimYear,ncol=MaxAge)
for(i in 1:SimYear)
 MovementS[i,]	<-MaxMovingS[i]/(1+exp(-1*log(19)*(Ages-Move50s[i])/(Move95s[i]-Move50s[i])))
}

if(sum(MovementS)>0 | sum(MovementN)>0) 
 {print("Movement is occuring between populations")}

#==recruitment parameters==
steepnessN	<-CleanInput(out$OM$steepnessN,SimYear)
cat(steepnessN,"\t # steepnessN \n",file=echofile,append=T)
sigmaRn	<-CleanInput(out$OM$sigmaRn,SimYear)
cat(sigmaRn,"\t # sigmaRn \n",file=echofile,append=T)
RzeroN	<-CleanInput(out$OM$RzeroN,SimYear)
cat(RzeroN,"\t # RzeroN \n",file=echofile,append=T)

steepnessS	<-steepnessN
sigmaRs	<-sigmaRn
RzeroS	<-RzeroN

if(TwoPop>0)
{
steepnessS	<-CleanInput(out$OM$steepnessS,SimYear)
sigmaRs	<-CleanInput(out$OM$sigmaRs,SimYear)
RzeroS	<-CleanInput(out$OM$RzeroS,SimYear)
}

RecErrN	<-matrix(rnorm(SimYear*Nsim,0,sigmaRn),ncol=SimYear,byrow=T)
cat("# RecErrN \n",file=echofile,append=T)
write.table(row.names=F,col.names=F,RecErrN,file=echofile,append=T)
RecErrS	<-matrix(rnorm(SimYear*Nsim,0,sigmaRs),ncol=SimYear,byrow=T)

#==historic fishing mortality
HistoricalFn	<-CleanInput(out$OM$HistoricalFn,SimYear)[1:InitYear]
cat(HistoricalFn,"\t # HistoricalFn \n",file=echofile,append=T)
PastFsdN		<-out$OM$PastFsdN
HarvestControlN	<-out$OM$HarvestControlN
ConstantCatchN	<-out$OM$ConstantCatchN
ConstantFn		<-out$OM$ConstantFn
HCalphaN		<-out$OM$HCalphaN
HCbetaN		<-out$OM$HCbetaN

HistoricalFs	<-HistoricalFn
PastFsdS		<-PastFsdN
HarvestControlS	<-HarvestControlN
ConstantCatchS	<-ConstantCatchN
ConstantFs		<-ConstantFn
HCalphaS		<-HCalphaN
HCbetaS		<-HCbetaN

if(TwoPop>0)
{
HistoricalFs	<-CleanInput(out$OM$HistoricalFs,SimYear)[1:InitYear]
PastFsdS		<-out$OM$PastFsdS
HarvestControlS	<-out$OM$HarvestControlS
ConstantCatchS	<-out$OM$ConstantCatchS
ConstantFs		<-out$OM$ConstantFs
HCalphaS		<-out$OM$HCalphaS
HCbetaS		<-out$OM$HCbetaS
}

start_assessment<-out$OM$start_assessment
cat(start_assessment,"\t # start_assessment \n",file=echofile,append=T)
SmallNum		<-out$OM$SmalNum
InitSmooth		<-out$OM$InitSmooth	
FmortPen		<-out$OM$FmortPen
RecruitPen		<-out$OM$RecruitPen	
Mpenalty		<-out$OM$Mpenalty	
Growthpenalty	<-out$OM$Growthpenalty	
SelPenalty		<-out$OM$SelPenalty

EstM			<-out$OM$EstM
TimeVaryM		<-out$OM$TimeVaryM	
EstGrowthK		<-out$OM$EstGrowthK
TimeVaryGrowthK	<-out$OM$TimeVaryGrowthK
EstLinf		<-out$OM$EstLinf
TimeVaryLinf	<-out$OM$TimeVaryLinf
TimeVarySel50	<-out$OM$TimeVarySel50
TimeVarySel95	<-out$OM$TimeVarySel95
ProjectTimeVary	<-out$OM$ProjectTimeVary
InitValSig		<-out$OM$InitValSig

InitBzeroMod	<-out$OM$InitBzeroMod
InitGrowthRate	<-out$OM$InitGrowthRate
estInit       <-out$OM$estInit
InitBioProd   <-out$OM$InitBioProd

#===================================================
#==Virgin numbers at age, biomass, recruitment
#===================================================
cat("\n\n # Initialised Population \n",file=echofile,append=T)
VirInitN		<-initialN(Rzero=RzeroN[1],NatM=NatMn[1],inAge=MaxAge)
cat("# VirInitN \n",file=echofile,append=T)
write.table(row.names=F,col.names=F,VirInitN,file=echofile,append=T)
VirInitS		<-initialN(Rzero=RzeroS[1],NatM=NatMs[1],inAge=MaxAge)

VirBioN		<-sum(VirInitN*matureN[1,]*WeightAtAgeN[1,])
cat("\n # VirBioN \n",file=echofile,append=T)
write.table(row.names=F,col.names=F,VirBioN,file=echofile,append=T)
VirBioS		<-sum(VirInitS*matureS[1,]*WeightAtAgeS[1,])

ExploitBioN		<-sum(VirInitN*vulnN[1,]*WeightAtAgeN[1,])
cat("\n # ExploitBioN \n",file=echofile,append=T)
write.table(row.names=F,col.names=F,ExploitBioN,file=echofile,append=T)
ExploitBioS		<-sum(VirInitS*vulnS[1,]*WeightAtAgeS[1,])

png(file.path(PlotFolder,paste0("LifeHistory_",CreateFolderName,".png")),res=1200,units='in',height=7.5,width=7.5)
par(mfrow=c(4,4),mar=c(3,3,0,0),oma=c(1,3,1,1))
PlotLifeHistory(LenAtAgeN,LenAtAgeS,matureN,matureS,vulnN,vulnS,survSelN,survSelS,WeightAtAgeN,
	WeightAtAgeS,MovementN,MovementS,NatMs,NatMn,VirBioN,VirBioS,RzeroN,RecErrN,steepnessN,steepnessS,
	RzeroS,RecErrS,sigmaRn,sigmaRs,HistoricalFn,HistoricalFs,SimYear,MaxAge)
dev.off()
cat("\n ####### PLOTTED LIFE HISTORY ####### \n",file=echofile,append=T)

#=========================================================================
# INITALIZE THE POPULATION
#==========================================================================
tempNn		<-array(dim=c(InitYear,MaxAge,Nsim))
tempNn[1,,]		<-VirInitN
tempNs		<-array(dim=c(InitYear,MaxAge,Nsim))
tempNs[1,,]		<-VirInitS
tempCatchN		<-matrix(ncol=InitYear,nrow=Nsim)
tempCatchS		<-matrix(ncol=InitYear,nrow=Nsim)
tempRecN		<-matrix(ncol=InitYear,nrow=Nsim)
tempRecS		<-matrix(ncol=InitYear,nrow=Nsim)
tempCatchAtAgeN	<-array(dim=c(InitYear,MaxAge,Nsim))
tempCatchAtAgeS	<-array(dim=c(InitYear,MaxAge,Nsim))

#==Make a matrix of past Fs based on a sd and an input time series==
HistoricalFsInit	<-matrix(HistoricalFs,ncol=InitYear,nrow=Nsim,byrow=T)
HistoricalFnInit	<-matrix(HistoricalFn,ncol=InitYear,nrow=Nsim,byrow=T)
FerrN			<-matrix(rnorm(InitYear*Nsim,1,PastFsdN),ncol=InitYear)
#cat("\n",file=echofile,append=T)
cat("\n # FerrN \n",file=echofile,append=T)
write.table(row.names=F,col.names=F,FerrN,file=echofile,append=T)
FerrS			<-matrix(rnorm(InitYear*Nsim,1,PastFsdS),ncol=InitYear)
HistoricalFsIn  <-HistoricalFsInit*FerrN
#cat("\n",file=echofile,append=T)
cat("\n # HistoricalFsIn \n",file=echofile,append=T)
write.table(row.names=F,col.names=F,HistoricalFsIn,file=echofile,append=T)
HistoricalFnIn  <-HistoricalFnInit*FerrS

for(k in 1:Nsim)
{
 for (j in 2:InitYear)
 {
  for (i in 2:(MaxAge-1))
  {
   tempNn[j,i,k]		<-tempNn[j-1,i-1,k]*exp(-HistoricalFnIn[k,j]*vulnN[1,i-1])*exp(-NatMn[j])
   tempNs[j,i,k]		<-tempNs[j-1,i-1,k]*exp(-HistoricalFsIn[k,j]*vulnS[1,i-1])*exp(-NatMs[j])
  }
   tempNn[j,MaxAge,k]	<-(tempNn[j-1,(MaxAge-1),k])*exp(-HistoricalFnIn[k,j]*vulnN[j,MaxAge])*exp(-NatMn[j])+ tempNn[j-1,MaxAge,k]*exp(-HistoricalFnIn[k,j]*vulnN[j,MaxAge])*exp(-NatMn[j])
   tempNs[j,MaxAge,k]	<-(tempNs[j-1,(MaxAge-1),k])*exp(-HistoricalFsIn[k,j]*vulnS[j,MaxAge])*exp(-NatMs[j])+ tempNs[j-1,MaxAge,k]*exp(-HistoricalFsIn[k,j]*vulnS[j,MaxAge])*exp(-NatMs[j])

   moveFromN		<-tempNn[j,,k]*MovementN[j,]
   moveFromS		<-tempNs[j,,k]*MovementS[j,]

   tempNn[j,,k]		<-tempNn[j,,k]-moveFromN+moveFromS
   tempNs[j,,k]		<-tempNs[j,,k]-moveFromS+moveFromN

   EggsN			<-sum(tempNn[j-1,,k]*matureN[j,]*WeightAtAgeN[j,])
   EggsS			<-sum(tempNs[j-1,,k]*matureS[j,]*WeightAtAgeS[j,])

   tempNn[j,1,k]		<-Recruitment(EggsIN=EggsN,steepnessIN=steepnessN[j],RzeroIN=RzeroN[j],RecErrIN=RecErrN[k,j],recType="BH",NatMin=NatMn[j],
							vulnIN=vulnN[j,],matureIN=matureN[j,],weightIN=WeightAtAgeN[j,],LenAtAgeIN=LenAtAgeN[j,],MaxAge=MaxAge,sigmaRin=sigmaRn[j])
   tempNs[j,1,k]		<-Recruitment(EggsIN=EggsS,steepnessIN=steepnessS[j],RzeroIN=RzeroS[j],RecErrIN=RecErrS[k,j],recType="BH",NatMin=NatMs[j],
							vulnIN=vulnS[j,],matureIN=matureS[j,],weightIN=WeightAtAgeS[j,],LenAtAgeIN=LenAtAgeS[j,],MaxAge=MaxAge,sigmaRin=sigmaRs[j])
   tempRecN[k,j]		<-tempNn[j,1,k]
   tempRecS[k,j]		<-tempNs[j,1,k]

   tempCatchAtAgeN[j,,k]<-((vulnN[j,]*HistoricalFnIn[k,j])/(vulnN[j,]*HistoricalFnIn[k,j]+NatMn[j])) * (1-exp(-(vulnN[j,]*HistoricalFnIn[k,j]+NatMn[j]))) * tempNn[j-1,,k]
   tempCatchN[k,j]	<-sum(tempCatchAtAgeN[j,,k]*WeightAtAgeN[j,])

   tempCatchAtAgeS[j,,k]<-((vulnS[j,]*HistoricalFsIn[k,j])/(vulnS[j,]*HistoricalFsIn[k,j]+NatMs[j])) * (1-exp(-(vulnS[j,]*HistoricalFsIn[k,j]+NatMs[j]))) * tempNs[j-1,,k]
   tempCatchS[k,j]	<-sum(tempCatchAtAgeS[j,,k]*WeightAtAgeS[j,])
 }
}
#===============================================================
# BEGIN SIMULATION OF ASSESSMENT AND HARVEST
#===============================================================
#==tempNn and tempNs from above are the starting points
cat("\n\n",file=echofile,append=T)
cat("## Full Population  \n",file=echofile,append=T)
projNn	<-array(dim=c(SimYear,MaxAge,Nsim))
projNs	<-array(dim=c(SimYear,MaxAge,Nsim))

projCatchAtAgeN	<-array(dim=c(SimYear,MaxAge,Nsim))
projCatchAtAgeS	<-array(dim=c(SimYear,MaxAge,Nsim))

projCatchN	<-matrix(nrow=Nsim,ncol=SimYear)
projCatchS	<-matrix(nrow=Nsim,ncol=SimYear)
projSSBn	<-matrix(nrow=Nsim,ncol=SimYear)
projSSBs	<-matrix(nrow=Nsim,ncol=SimYear)
projExpBn	<-matrix(nrow=Nsim,ncol=SimYear)
projExpBs	<-matrix(nrow=Nsim,ncol=SimYear)
projSurvN	<-matrix(nrow=Nsim,ncol=SimYear)
projSurvS	<-matrix(nrow=Nsim,ncol=SimYear)
projRecN	<-matrix(nrow=Nsim,ncol=SimYear)
projRecS	<-matrix(nrow=Nsim,ncol=SimYear)
projFmortN	<-matrix(nrow=Nsim,ncol=SimYear)
projFmortS	<-matrix(nrow=Nsim,ncol=SimYear)

projCatLenFreqN	<-array(dim=c(SimYear,LengthBinN,Nsim))
projCatLenFreqS	<-array(dim=c(SimYear,LengthBinN,Nsim))
projSurvLenFreqN	<-array(dim=c(SimYear,LengthBinN,Nsim))
projSurvLenFreqS	<-array(dim=c(SimYear,LengthBinN,Nsim))

#=storage for the true quantities
trueRecN		<-matrix(ncol=SimYear,nrow=Nsim)
trueCatchN		<-matrix(ncol=SimYear,nrow=Nsim)

trueRecS		<-matrix(ncol=SimYear,nrow=Nsim)
trueCatchS		<-matrix(ncol=SimYear,nrow=Nsim)

trueFmortN		<-matrix(ncol=SimYear,nrow=Nsim)
for(x in 1:Nsim)
 trueFmortN[x,1:InitYear]<-HistoricalFn
trueFmortS		<-matrix(ncol=SimYear,nrow=Nsim)
for(x in 1:Nsim)
 trueFmortS[x,1:InitYear]<-HistoricalFn

trueSpbioN		<-matrix(ncol=SimYear,nrow=Nsim)
trueSurvIndN	<-matrix(ncol=SimYear,nrow=Nsim)
trueCPUEindN	<-matrix(ncol=SimYear,nrow=Nsim)

trueSpbioS		<-matrix(ncol=SimYear,nrow=Nsim)
trueSurvIndS	<-matrix(ncol=SimYear,nrow=Nsim)
trueCPUEindS	<-matrix(ncol=SimYear,nrow=Nsim)

trueF35   		<-rep(0,SimYear)
trueSBPR35 		<-rep(0,SimYear)
trueB35		    <-matrix(ncol=SimYear,nrow=Nsim)
#========================================================================
# FIND REFERENCE POINTS FOR THE POPULATION 
#========================================================================
if(PlotYieldCurve==1)
{
pdf(file.path(PlotFolder,paste0("Plot_init_YieldCurve_",CreateFolderName,".pdf")))
SearchFmort		<-seq(0.01,3*NatMn[1],(NatMn[1]-0.01)/100)
SearchYield		<-rep(0,length(SearchFmort))
SearchBiomass	<-rep(0,length(SearchFmort))
print("Plotting yield curve")
for(p in 1:length(SearchFmort))
{
tempOut<-ProjPopDym(fmortN=SearchFmort[p],MaxAge=MaxAge,vulnN=vulnN,
                     	vulnS=vulnS,NatMn=NatMn,NatMs=NatMs,matureN=matureN,matureS=matureS,
				WeightAtAgeN=WeightAtAgeN, WeightAtAgeS=WeightAtAgeS,steepnessN=steepnessN,
				steepnessS=steepnessS,RzeroN=RzeroN,RzeroS=RzeroS,LenAtAgeN=LenAtAgeN,
				LenAtAgeS=LenAtAgeS,sigmaRn=sigmaRn,sigmaRs=sigmaRs,RefYear=SimYear,
        MovementN=MovementN,MovementS=MovementS)

SearchYield[p]<-tempOut[[1]]
SearchBiomass[p]<-tempOut[[2]]
}
print("Yield curve plotted")
par(mfrow=c(1,2))
plot(SearchYield~SearchFmort)
plot(SearchYield~SearchBiomass)
dev.off()
}

#==============================================================
# calculate the catch proportion at (length) for assessment
#==============================================================
 tempCatAtLenN<-array(dim=c(ncol(LenAtAgeN),LengthBinN,Nsim))
 tempCatAtLenS<-array(dim=c(ncol(LenAtAgeN),LengthBinN,Nsim))

for(x in 1:Nsim)
 for(y in 2:InitYear)
 {
 #==make length frequencies for catch==
 for(w in 1:ncol(LenAtAgeN))
 {
  probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeN[y,w],sd=GrowthSDn[y])
  ProbN<-probtemp/sum(probtemp)
  tempCatAtLenN[w,,x]<-tempCatchAtAgeN[y,w,x]*ProbN
  probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeS[y,w],sd=GrowthSDs[y])
  ProbS<-probtemp/sum(probtemp)
  tempCatAtLenS[w,,x]<-tempCatchAtAgeS[y,w,x]*ProbS
 }
 projCatLenFreqN[y,,x]<-apply(tempCatAtLenN[,,x],2,sum)
 projCatLenFreqS[y,,x]<-apply(tempCatAtLenS[,,x],2,sum)
}
cat("\n",file=echofile,append=T)
cat(paste0("# projCatLenFreqN #array! \n"),file=echofile,append=T)
for(n in 1:Nsim) {
  cat(paste0("# projCatLenFreqN sim: ",n," \n"),file=echofile,append=T)
  write.table(row.names=F,col.names=F,projCatLenFreqN[,,n],file=echofile,append=T)
}

if(AssessmentData==1)
 {
 pdf(file.path(PlotFolder,paste0("CatchLengthProportions_",CreateFolderName,".pdf")))
 par(mfrow=c(ceiling(sqrt(InitYear)),ceiling(sqrt(InitYear))),mar=c(.1,.1,.1,.1))
 for(i in 2:InitYear)
  barplot(rbind(projCatLenFreqN[i,,1],projCatLenFreqS[i,,1]),beside=T,xaxt='n',yaxt='n')
 dev.off()
 }

#==============================================================
# calculate the survey proportion at length for assessment
#==============================================================
 tempSurvAtLenN<-array(dim=c(ncol(LenAtAgeN),LengthBinN,Nsim))
 tempSurvAtLenS<-array(dim=c(ncol(LenAtAgeN),LengthBinN,Nsim))

for(x in 1:Nsim)
 for(y in 2:InitYear)
 {
 #==make length frequencies for s==
 for(w in 1:ncol(LenAtAgeS))
 {
  probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeN[y,w],sd=GrowthSDn[y])
  ProbN<-probtemp/sum(probtemp)
  tempSurvAtLenN[w,,x]<-tempNn[y,w,x]*survSelN[y,w]*ProbN
  probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeS[y,w],sd=GrowthSDs[y])
  ProbS<-probtemp/sum(probtemp)
  tempSurvAtLenS[w,,x]<-tempNs[y,w,x]*survSelS[y,w]*ProbS
 }
 projSurvLenFreqN[y,,x]<-apply(tempSurvAtLenN[,,x],2,sum)
 projSurvLenFreqS[y,,x]<-apply(tempSurvAtLenS[,,x],2,sum)
}
cat("\n",file=echofile,append=T)
cat(paste0("# projSurvLenFreqN #array!"," \n"),file=echofile,append=T)
for(n in 1:Nsim) {
  cat(paste0("# projSurvLenFreqN sim:",n," \n"),file=echofile,append=T)  
  write.table(row.names=F,col.names=F,projSurvLenFreqN[,,n],file=echofile,append=T)
}

if(AssessmentData==1)
 {
 pdf(file.path(PlotFolder,paste0("SurveyLengthProportions_",CreateFolderName,".pdf")))
 par(mfrow=c(ceiling(sqrt(InitYear)),ceiling(sqrt(InitYear))),mar=c(.1,.1,.1,.1))
 for(i in 2:InitYear)
  barplot(rbind(projSurvLenFreqN[i,,1],projSurvLenFreqS[i,,1]),beside=T,xaxt='n',yaxt='n')
dev.off()
 }

#==transfer historical time series from above
for(x in 1:Nsim)
{
projNn[1:InitYear,,x]			<-tempNn[,,x]
projNs[1:InitYear,,x]			<-tempNs[,,x]
projCatchN[x,1:InitYear]		<-tempCatchN[x,]
projCatchS[x,1:InitYear]		<-tempCatchS[x,]
projRecN[x,1:InitYear]			<-tempRecN[x,]
projRecS[x,1:InitYear]			<-tempRecS[x,]
projFmortN[x,1:InitYear]		<-HistoricalFnInit[x,]
projFmortS[x,1:InitYear]		<-HistoricalFsInit[x,]
projSurvN[x,1:InitYear]			<-apply(tempNn[,,x]*survSelN[1:InitYear,]*WeightAtAgeN[1:InitYear,],1,sum)
projSurvS[x,1:InitYear]			<-apply(tempNs[,,x]*survSelS[1:InitYear,]*WeightAtAgeS[1:InitYear,],1,sum)

for(y in 1:InitYear)
{
 projSSBn[x,y]	<-sum(projNn[y,,x]*matureN[y,]*WeightAtAgeN[y,])
 projSSBs[x,y]	<-sum(projNs[y,,x]*matureS[y,]*WeightAtAgeS[y,])
 projExpBn[x,y]	<-sum(projNn[y,,x]*vulnN[y,]*WeightAtAgeN[y,])
 projExpBs[x,y]	<-sum(projNs[y,,x]*vulnS[y,]*WeightAtAgeS[y,])
}
}
cat(paste0("# projNn #array!"," \n"),file=echofile,append=T)
for(n in 1:Nsim) {
  cat(paste0("# projNn sim:",n," \n"),file=echofile,append=T)
  write.table(row.names=F,col.names=F,projNn[,,n],file=echofile,append=T)
}
cat("\n",file=echofile,append=T)
cat("# projRecN \n",file=echofile,append=T)
write.table(row.names=F,col.names=F,projRecN,file=echofile,append=T)
cat("\n",file=echofile,append=T)
cat("# projFmortN \n",file=echofile,append=T)
write.table(row.names=F,col.names=F,projFmortN,file=echofile,append=T)
cat("\n",file=echofile,append=T)
cat("# projSSBn \n",file=echofile,append=T)
write.table(row.names=F,col.names=F,projSSBn,file=echofile,append=T)
cat("\n",file=echofile,append=T)
cat("# projExpBn \n",file=echofile,append=T)
write.table(row.names=F,col.names=F,projExpBn,file=echofile,append=T)
cat("\n",file=echofile,append=T)
cat("# projSurvN \n",file=echofile,append=T)
write.table(row.names=F,col.names=F,projSurvN,file=echofile,append=T)

projCatchN[,1]		<-0
cat("\n",file=echofile,append=T)
cat("# projCatchN \n",file=echofile,append=T)
write.table(row.names=F,col.names=F,projCatchN,file=echofile,append=T)

projCatchS[,1]		<-0

#==fill in true storage arrays
for(x in 1:Nsim)
{
trueRecN[x,1:InitYear]			<-tempRecN[x,]
trueCatchN[x,1:InitYear]		<-tempCatchN[x,]		
trueSpbioN[x,1:InitYear]		<-projSSBn[x,1:InitYear]
trueSurvIndN[x,1:InitYear]		<-projSurvN[x,1:InitYear]
trueCPUEindN[x,1:InitYear]		<-projExpBn[x,1:InitYear]
trueRecS[x,1:InitYear]			<-tempRecS[x,]	
trueCatchS[x,1:InitYear]		<-tempCatchS[x,]		
trueSpbioS[x,1:InitYear]		<-projSSBs[x,1:InitYear]	
trueSurvIndS[x,1:InitYear]		<-projSurvS[x,1:InitYear]
trueCPUEindS[x,1:InitYear]		<-projExpBs[x,1:InitYear]
}

CatchErrorN		<-matrix(rnorm(Nsim*SimYear,0,CatchCVn),nrow=Nsim,ncol=SimYear)
cat("\n # CatchErrorN \n",file=echofile,append=T)
write.table(row.names=F,col.names=F,CatchErrorN,file=echofile,append=T)
CatchErrorS		<-matrix(rnorm(Nsim*SimYear,0,CatchCVs),nrow=Nsim,ncol=SimYear)
CatchAssessN	<-projCatchN*exp(CatchErrorN-(CatchCVn^2/2))
cat("\n # CatchAssessN \n",file=echofile,append=T)
write.table(row.names=F,col.names=F,CatchAssessN,file=echofile,append=T)
CatchAssessS	<-projCatchS*exp(CatchErrorS-(CatchCVs^2/2))

CPUEErrorN		<-matrix(rnorm(Nsim*SimYear,0,IndexCVn),nrow=Nsim,ncol=SimYear)
cat("\n # CPUEErrorN \n",file=echofile,append=T)
write.table(row.names=F,col.names=F,CPUEErrorN,file=echofile,append=T)
CPUEErrorS		<-matrix(rnorm(Nsim*SimYear,0,IndexCVs),nrow=Nsim,ncol=SimYear)
CPUEAssessN		<-projExpBn*exp(CPUEErrorN-(IndexCVn^2/2))
cat("\n # CPUEAssessN \n",file=echofile,append=T)
write.table(row.names=F,col.names=F,CPUEAssessN,file=echofile,append=T)
CPUEAssessS		<-projExpBs*exp(CPUEErrorS-(IndexCVs^2/2))

SurvErrorN		<-matrix(rnorm(Nsim*SimYear,0,IndexCVn),nrow=Nsim,ncol=SimYear)
cat("\n # SurvErrorN \n",file=echofile,append=T)
write.table(row.names=F,col.names=F,SurvErrorN,file=echofile,append=T)
SurvErrorS		<-matrix(rnorm(Nsim*SimYear,0,IndexCVs),nrow=Nsim,ncol=SimYear)
SurvAssessN		<-projSurvN*exp(SurvErrorN-(IndexCVn^2/2))
cat("\n # SurveAssessN \n",file=echofile,append=T)
write.table(row.names=F,col.names=F,SurvAssessN,file=echofile,append=T)
SurvAssessS		<-projSurvS*exp(SurvErrorS-(IndexCVs^2/2))

cat(paste0("Input Ended: ",Sys.time(),"\n"),file=echofile,append=T)

#==storage for management and assessment quantities
FMSY		<-matrix(nrow=Nsim,ncol=SimYear)
BMSY		<-matrix(nrow=Nsim,ncol=SimYear)
TAC		<-matrix(nrow=Nsim,ncol=SimYear)
trueTAC	<-matrix(nrow=Nsim,ncol=SimYear)
trueOFL	<-matrix(nrow=Nsim,ncol=SimYear)
CurBio	<-matrix(nrow=Nsim,ncol=SimYear)
EstBio	<-matrix(nrow=Nsim,ncol=SimYear)
Converge	<-matrix(nrow=Nsim,ncol=SimYear)
est_init_B <-matrix(nrow=Nsim,ncol=SimYear)

#============================================
# calculate 'true' reference points
# figure out a way to determine when to do this for all sims and when for just 1
# as long as population parametesr are from the same control file, this should be the same 
# for all populations regardless of fishing history or recruitment
#===================================
FindFMSYin<-1
if(FindFMSYin==1)
{
  #==Find reference points for the initial year of projection
  x			<-0.3
  FmsyOut		<-nlminb(x,FindFMSY,MaxAge=MaxAge,VirInitN=VirInitN,VirInitS=VirInitS,vulnN=vulnN,
                    vulnS=vulnS,NatMn=NatMn,NatMs=NatMs,matureN=matureN,matureS=matureS,
                    WeightAtAgeN=WeightAtAgeN, WeightAtAgeS=WeightAtAgeS,steepnessN=steepnessN,
                    steepnessS=steepnessS,RzeroN=RzeroN,RzeroS=RzeroS,LenAtAgeN=LenAtAgeN,LenAtAgeS=LenAtAgeS,
                    sigmaRn=sigmaRn,sigmaRs=sigmaRs,RefYear=InitYear,MovementN=MovementN,MovementS=MovementS)
  trueFMSYbegin	<-FmsyOut$par
  trueUMSYbegin	<- 1-(exp(-trueFMSYbegin))
  trueBMSYbegin	<-ProjPopDym(fmortN=trueFMSYbegin,MaxAge=MaxAge,vulnN=vulnN,
                             vulnS=vulnS,NatMn=NatMn,NatMs=NatMs,matureN=matureN,matureS=matureS,
                             WeightAtAgeN=WeightAtAgeN, WeightAtAgeS=WeightAtAgeS,steepnessN=steepnessN,
                             steepnessS=steepnessS,RzeroN=RzeroN,RzeroS=RzeroS,LenAtAgeN=LenAtAgeN,LenAtAgeS=LenAtAgeS,
                             sigmaRn=sigmaRn,sigmaRs=sigmaRs,RefYear=InitYear,MovementN=MovementN,MovementS=MovementS)[[2]]
  trueMSYbegin	<- -FmsyOut$objective
  
  #==find the reference points associated with the final year of projection
  FmsyOut		<-nlminb(x,FindFMSY,MaxAge=MaxAge,VirInitN=VirInitN,VirInitS=VirInitS,vulnN=vulnN,
                    vulnS=vulnS,NatMn=NatMn,NatMs=NatMs,matureN=matureN,matureS=matureS,
                    WeightAtAgeN=WeightAtAgeN, WeightAtAgeS=WeightAtAgeS,steepnessN=steepnessN,
                    steepnessS=steepnessS,RzeroN=RzeroN,RzeroS=RzeroS,LenAtAgeN=LenAtAgeN,LenAtAgeS=LenAtAgeS,
                    sigmaRn=sigmaRn,sigmaRs=sigmaRs,RefYear=SimYear,MovementN=MovementN,MovementS=MovementS)
  trueFMSYend		<-FmsyOut$par
  trueUMSYend		<- 1-(exp(-trueFMSYend))
  trueBMSYend		<-ProjPopDym(fmortN=trueFMSYend,MaxAge=MaxAge,vulnN=vulnN,
                            vulnS=vulnS,NatMn=NatMn,NatMs=NatMs,matureN=matureN,matureS=matureS,
                            WeightAtAgeN=WeightAtAgeN, WeightAtAgeS=WeightAtAgeS,steepnessN=steepnessN,
                            steepnessS=steepnessS,RzeroN=RzeroN,RzeroS=RzeroS,LenAtAgeN=LenAtAgeN,
                            LenAtAgeS=LenAtAgeS,sigmaRn=sigmaRn,sigmaRs=sigmaRs,RefYear=SimYear,
                            MovementN=MovementN,MovementS=MovementS)[[2]]
  trueMSYend		<- -FmsyOut$objective
}

#======================================
# calculate reference point proxies
# MIGHT NEED TWO OF THESE
#======================================
ConstRec	<-100000
outsF35in	<-FindF35(MaxAge=MaxAge,vulnN=vulnN,
                    vulnS=vulnS,NatMn=NatMn,NatMs=NatMs,matureN=matureN,matureS=matureS,
                    WeightAtAgeN=WeightAtAgeN, WeightAtAgeS=WeightAtAgeS,steepnessN=steepnessN,
                    steepnessS=steepnessS,RzeroN=RzeroN,RzeroS=RzeroS,LenAtAgeN=LenAtAgeN,
                    LenAtAgeS=LenAtAgeS,inRec=ConstRec,RefYear=InitYear,MovementN=MovementN,MovementS=MovementS)

trueSBPR35[InitYear]<-outsF35in[[2]]
trueF35[InitYear]	  <-outsF35in[[1]]
trueB35[,InitYear]  <-trueSBPR35[InitYear]*apply(trueRecN[,start_assessment:SimYear],1,mean,na.rm=T)

#=================================================================================
#==BEGIN PROJECTIONS===========================================================
#=================================================================================

for(z in 1:Nsim)
{
 for(y in (InitYear+1):SimYear)
 {
  if(z==1)
  {
    # this calculates the F35 and SBPR35
    # make it so this only does it if the things that change F35 and SBPR35 change
    # this can be done with an if statement that calls a function that checksif the values are all the same (or columns are all the same)
    outsF35in	<-FindF35(MaxAge=MaxAge,vulnN=vulnN,
                        vulnS=vulnS,NatMn=NatMn,NatMs=NatMs,matureN=matureN,matureS=matureS,
                        WeightAtAgeN=WeightAtAgeN, WeightAtAgeS=WeightAtAgeS,steepnessN=steepnessN,
                        steepnessS=steepnessS,RzeroN=RzeroN,RzeroS=RzeroS,LenAtAgeN=LenAtAgeN,
                        LenAtAgeS=LenAtAgeS,inRec=ConstRec,RefYear=y,MovementN=MovementN,MovementS=MovementS)
    
    trueSBPR35[y]<-outsF35in[[2]]
    trueF35[y]	 <-outsF35in[[1]]
  }
 #==pull data for assessment
 CatchDataN	<- CatchAssessN[z,!is.na(CatchAssessN[z,])]
 CPUEDataN	<- CPUEAssessN[z,!is.na(CPUEAssessN[z,])]
 SurvDataN	<- SurvAssessN[z,!is.na(SurvAssessN[z,])] 

 CatchDataS	<- CatchAssessS[z,!is.na(CatchAssessS[z,])]
 CPUEDataS	<- CPUEAssessS[z,!is.na(CPUEAssessS[z,])] 
 SurvDataS	<- SurvAssessS[z,!is.na(SurvAssessS[z,])] 

 #=================================================================================
 #===ASSESS THE STOCK, CHOOSE A METHOD IN CTL FILE=================================
 #=================================================================================

 #==Production model==
 if(AssessmentType == 1)
 {
 setwd(MSEdir)
 if(!dir.exists(CreateFolderName)) dir.create(CreateFolderName)
 setwd(CreateFolderName)
 if(y==(InitYear+1))
 {
  x		<-c(log(InitBzeroMod*(VirBioN)),InitGrowthRate)	
  if(estInit==1)
    x		<-c(log(InitBzeroMod*(VirBioN)),InitGrowthRate,InitBioProd)	
 }
  if(y>(InitYear+1))
  x		<-outs$par

 inCatch	<-CatchDataN[start_assessment:(y-1)]
 inCPUE	<-CPUEDataN[start_assessment:(y-1)]
 outs		<-nlminb(start=x,objective=ProdMod,CatchData=inCatch,IndexData=inCPUE,estInit=estInit)
# outs <- optim(par=x,fn=ProdMod,CatchData=inCatch,IndexData=inCPUE,estInit=estInit)
 Converge[z,y]<-outs$convergence
 PredBio	<-ProdModPlot(outs$par,inCatch,inCPUE,plots=EstimationPlots,estInit=estInit)
 FMSY[z,y] 	<-outs$par[2]/2
 BMSY[z,y]	<-exp(outs$par[1])/2
 if(estInit==1)
  est_init_B[z,y]<-outs$par[3]
 
 MSY		<-BMSY[z,y]*FMSY[z,y]
 CurBio[z,y]<-PredBio[length(PredBio)-1]
 inBio	<-projExpBn[z,y-1]
 TAC[z,y]	<-HarvestControlRule(FMSY=FMSY[z,y],BMSY=BMSY[z,y],ExploitBio=CurBio[z,y],SpawnBio=CurBio[z,y],
			alpha=HCalphaN,beta=HCbetaN,HarvestControl=HarvestControlN,ConstantCatch=ConstantCatchN,
  			ConstantF=ConstantFn)
 trueTAC[z,y]<-HarvestControlRule(FMSY=trueFMSYend,BMSY=trueBMSYend,ExploitBio=inBio,SpawnBio=inBio,
			alpha=HCalphaN,beta=HCbetaN,HarvestControl=HarvestControlN,ConstantCatch=ConstantCatchN,
  			ConstantF=ConstantFn)

 if(EstimationPlots==1 & AssessmentType == 1)
 {
  lines(BMSY[z,],col=3,lty=3)
  legend("topright",bty='n',col=c(NA,1,2,3),pch=c(15,NA,NA,NA),lty=c(NA,1,2,3),
          legend=c("Index obs","Index est","Catch","BMSY"))
 }

 if(y==SimYear)
 {
  EstBio[z,start_assessment:ncol(EstBio)]<-PredBio
 }
 if(y==SimYear & z==Nsim)
 {
 inFile<-paste("ProdOutputs.csv",sep="")
 file.create(inFile)

 cat("#True quantities","\n",file=inFile)
 cat("#","\n",file=inFile,append=TRUE)
 cat("#true Catch","\n",file=inFile,append=TRUE)
 for(m in 1:Nsim)
 cat(trueCatchN[m,],"\n",file=inFile,append=TRUE)
 cat("#true fishing mortality","\n",file=inFile,append=TRUE)
 for(m in 1:Nsim)
 cat(trueFmortN[m,],"\n",file=inFile,append=TRUE)
 cat("#true spawning biomass","\n",file=inFile,append=TRUE)
 for(m in 1:Nsim)
 cat(trueSpbioN[m,],"\n",file=inFile,append=TRUE)
 cat("#true CPUE","\n",file=inFile,append=TRUE)
 for(m in 1:Nsim)
 cat(trueCPUEindN[m,],"\n",file=inFile,append=TRUE)
 cat("#true total allowable catch","\n",file=inFile,append=TRUE)
 for(m in 1:Nsim)
 cat(trueTAC[m,],"\n",file=inFile,append=TRUE)
 cat("#true BMSY begin","\n",file=inFile,append=TRUE)
 cat(trueBMSYbegin,"\n",file=inFile,append=TRUE)
 cat("#true BMSY end","\n",file=inFile,append=TRUE)
 cat(trueBMSYend,"\n",file=inFile,append=TRUE)
 cat("#true FMSY begin","\n",file=inFile,append=TRUE)
 cat(trueFMSYbegin,"\n",file=inFile,append=TRUE)
 cat("#true FMSY end","\n",file=inFile,append=TRUE)
 cat(trueFMSYend,"\n",file=inFile,append=TRUE)
 cat("#est cpue ind","\n",file=inFile,append=TRUE)
 for(m in 1:Nsim)
 cat(EstBio[m,],"\n",file=inFile,append=TRUE)
 cat("#est spawning biomass by year (terminal year)","\n",file=inFile,append=TRUE)
 for(m in 1:Nsim)
 cat(CurBio[m,],"\n",file=inFile,append=TRUE)
 cat("#est total allowable catch","\n",file=inFile,append=TRUE)
 for(m in 1:Nsim)
 cat(TAC[m,],"\n",file=inFile,append=TRUE)
 cat("#est BMSY","\n",file=inFile,append=TRUE)
 for(m in 1:Nsim)
 cat(BMSY[m,],"\n",file=inFile,append=TRUE)
 cat("#est FMSY","\n",file=inFile,append=TRUE)
 for(m in 1:Nsim)
 cat(FMSY[m,],"\n",file=inFile,append=TRUE)
 if(estInit==1)
 {
   cat("#est initial biomass","\n",file=inFile,append=TRUE)
   for(m in 1:Nsim)
     cat(est_init_B[m,],"\n",file=inFile,append=TRUE) 
 }
 
 cat("#data cpue","\n",file=inFile,append=TRUE)
 for(m in 1:Nsim)
 cat(CPUEDataN,"\n",file=inFile,append=TRUE)
 }
}
#========================================
#==Projection based model--no assessment
#======================================= 
 if(AssessmentType == 0)
 {
 #==make folder if it isn't there
 setwd(MSEdir)
 if(!dir.exists(CreateFolderName)) dir.create(CreateFolderName)
 dir.create(file.path(CreateFolderName,z))
 IndSimFolder<-file.path(CreateFolderName,z,y)
 dir.create(IndSimFolder)

 #==write the true values
 setwd(IndSimFolder)
 file.create("TrueQuantities.DAT")
 cat("#True quantities","\n",file="TrueQuantities.DAT")
 cat("#","\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#recruitment","\n",file="TrueQuantities.DAT",append=TRUE)
 for(q in 1:Nsim)
 cat(trueRecN[q,],"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#Catch","\n",file="TrueQuantities.DAT",append=TRUE)
 for(q in 1:Nsim)
 cat(trueCatchN[q,],"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#fishing mortality","\n",file="TrueQuantities.DAT",append=TRUE)
 for(q in 1:Nsim)
 cat(trueFmortN[q,],"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#spawning biomass","\n",file="TrueQuantities.DAT",append=TRUE)
 for(q in 1:Nsim)
 cat(trueSpbioN[q,],"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#survey index","\n",file="TrueQuantities.DAT",append=TRUE)
 for(q in 1:Nsim)
 cat(trueSurvIndN[q,],"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#cpue index","\n",file="TrueQuantities.DAT",append=TRUE)
 for(q in 1:Nsim)
 cat(trueCPUEindN[q,],"\n",file="TrueQuantities.DAT",append=TRUE)
 #cat("#total allowable catch","\n",file="TrueQuantities.DAT",append=TRUE)
 #cat(trueTACn[q,],"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#","\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#BMSY begin","\n",file="TrueQuantities.DAT",append=TRUE)
 cat(trueBMSYbegin,"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#FMSY begin","\n",file="TrueQuantities.DAT",append=TRUE)
 cat(trueFMSYbegin,"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#BMSY end","\n",file="TrueQuantities.DAT",append=TRUE)
 cat(trueBMSYend,"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#FMSY end","\n",file="TrueQuantities.DAT",append=TRUE)
 cat(trueFMSYend,"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#B35","\n",file="TrueQuantities.DAT",append=TRUE)
 for(q in 1:Nsim)
 cat(trueB35[q,],"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#F35","\n",file="TrueQuantities.DAT",append=TRUE)
 cat(trueF35,"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#OFL","\n",file="TrueQuantities.DAT",append=TRUE)
 for(q in 1:Nsim)
 cat(trueOFL[q,],"\n",file="TrueQuantities.DAT",append=TRUE)

 #==calculate the TAC based on a SRR or a proxy, depending on data
 #==allow this to be a user defined function of any of the things that can be input
 TAC[z,y]	<-trueSpbio[z,y-1]*(1-exp(-ConstantF))
 }

 #======================================================
 #==Age-structured assessment model for setting the TAC
 #======================================================
 if(AssessmentType == 2)
 {
 #==make folder if it isn't there
 setwd(MSEdir)
 if(!dir.exists(CreateFolderName)) dir.create(CreateFolderName)
 dir.create(file.path(CreateFolderName,z))
 IndSimFolder<-file.path(CreateFolderName,z,y)
 dir.create(IndSimFolder)

 #==copy the .exe into it (where does this come from? github?)
 #==Better way of doing this??
 execfile <- file.path(GeMSdir,"GenAss",SimAssExec)
 if(!file.exists(GeMSdir)) stop("Executable file not found. Please set GeMSdir to GeMS directory. \n e.g. C:/GeMS")
 file.copy(from=execfile,to=IndSimFolder)

 #==write the true values
 #==Probably don't need this if it is stored in the MSE object
 setwd(IndSimFolder)
 file.create("TrueQuantities.DAT")
 cat("#True quantities","\n",file="TrueQuantities.DAT")
 cat("#","\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#recruitment","\n",file="TrueQuantities.DAT",append=TRUE)
 for(q in 1:Nsim)
 cat(trueRecN[q,],"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#Catch","\n",file="TrueQuantities.DAT",append=TRUE)
 for(q in 1:Nsim)
 cat(trueCatchN[q,],"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#fishing mortality","\n",file="TrueQuantities.DAT",append=TRUE)
 for(q in 1:Nsim)
 cat(trueFmortN[q,],"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#spawning biomass","\n",file="TrueQuantities.DAT",append=TRUE)
 for(q in 1:Nsim)
 cat(trueSpbioN[q,],"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#survey index","\n",file="TrueQuantities.DAT",append=TRUE)
 for(q in 1:Nsim)
 cat(trueSurvIndN[q,],"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#cpue index","\n",file="TrueQuantities.DAT",append=TRUE)
 for(q in 1:Nsim)
 cat(trueCPUEindN[q,],"\n",file="TrueQuantities.DAT",append=TRUE)
 #cat("#total allowable catch","\n",file="TrueQuantities.DAT",append=TRUE)
 #cat(trueTAC[q,],"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#","\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#BMSY begin","\n",file="TrueQuantities.DAT",append=TRUE)
 cat(trueBMSYbegin,"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#FMSY begin","\n",file="TrueQuantities.DAT",append=TRUE)
 cat(trueFMSYbegin,"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#BMSY end","\n",file="TrueQuantities.DAT",append=TRUE)
 cat(trueBMSYend,"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#FMSY end","\n",file="TrueQuantities.DAT",append=TRUE)
 cat(trueFMSYend,"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#B35in","\n",file="TrueQuantities.DAT",append=TRUE)
 for(q in 1:Nsim)
 cat(trueB35[q,],"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#F35in","\n",file="TrueQuantities.DAT",append=TRUE)
 cat(trueF35,"\n",file="TrueQuantities.DAT",append=TRUE)
 cat("#OFL","\n",file="TrueQuantities.DAT",append=TRUE)
 for(q in 1:Nsim)
 cat(trueOFL[q,],"\n",file="TrueQuantities.DAT",append=TRUE)

 #selAtAgeFunc(sel50n[1],VonKn[1],LinfN[1],t0n[1])
 #selAtAgeFunc(sel95n[1],VonKn[1],LinfN[1],t0n[1])
 #selAtAgeFunc(surv50n[1],VonKn[1],LinfN[1],t0n[1])
 #selAtAgeFunc(surv95n[1],VonKn[1],LinfN[1],t0n[1])

 #==.CTL file  !!!!!!!!!!!!!!!!!how should we deal with space?!!!!!!!!!!!!!!!!!!!!!!
 file.create("SimAss.CTL")
 cat("#Simulated assessment control file","\n",file="SimAss.CTL")
 cat("#","\n",file="SimAss.CTL",append=TRUE)
 cat("#weight parameters","\n",file="SimAss.CTL",append=TRUE)
 cat(c(alphaN[1],betaN[1]) ,"\n",file="SimAss.CTL",append=TRUE)
 cat("#how many years to average over timevarying process for projectsion","\n",file="SimAss.CTL",append=TRUE)
 cat(ProjectTimeVary,"\n",file="SimAss.CTL",append=TRUE)
 cat("#Estimate growth? if positive, indicates phase, no estimate if negative","\n",file="SimAss.CTL",append=TRUE)
 cat(EstGrowthK ,"\n",file="SimAss.CTL",append=TRUE)
 cat("#growth time-varying?","\n",file="SimAss.CTL",append=TRUE)
 cat(TimeVaryGrowthK ,"\n",file="SimAss.CTL",append=TRUE)
 cat("#Estimate growth? if positive, indicates phase, no estimate if negative","\n",file="SimAss.CTL",append=TRUE)
 cat(EstLinf ,"\n",file="SimAss.CTL",append=TRUE)
 cat("#growth time-varying?","\n",file="SimAss.CTL",append=TRUE)
 cat(TimeVaryLinf ,"\n",file="SimAss.CTL",append=TRUE)
 cat("#length parameters","\n",file="SimAss.CTL",append=TRUE)
 cat(c(VonKn[1],LinfN[1],t0n[1]) ,"\n",file="SimAss.CTL",append=TRUE)
 cat("#estimate natural mortality? if positive, number indicates phase, if negative, not estimated","\n",file="SimAss.CTL",append=TRUE)
 cat(EstM,"\n",file="SimAss.CTL",append=TRUE)
 cat("#natural mortality (used as starting value if estimated)","\n",file="SimAss.CTL",append=TRUE)
 cat(NatMn[length(NatMn)],"\n",file="SimAss.CTL",append=TRUE)
 cat("#Allow M to vary over time? 1 = yes, 0 = no","\n",file="SimAss.CTL",append=TRUE)
 cat(TimeVaryM,"\n",file="SimAss.CTL",append=TRUE) 
 cat("#selectivity 50% time-varying?","\n",file="SimAss.CTL",append=TRUE)
 cat(TimeVarySel50 ,"\n",file="SimAss.CTL",append=TRUE)
 cat("#selectivity 95% time-varying?","\n",file="SimAss.CTL",append=TRUE)
 cat(TimeVarySel95 ,"\n",file="SimAss.CTL",append=TRUE)
 cat("#maturity","\n",file="SimAss.CTL",append=TRUE)
 cat(c(mat50n[length(mat50n)],mat95n[length(mat95n)]) ,"\n",file="SimAss.CTL",append=TRUE)
 cat("#steepness","\n",file="SimAss.CTL",append=TRUE)
 cat(steepnessN[length(steepnessN)],"\n",file="SimAss.CTL",append=TRUE)
 cat("#small number","\n",file="SimAss.CTL",append=TRUE)
 cat(SmallNum,"\n",file="SimAss.CTL",append=TRUE)
 cat("#initial smoothness weight","\n",file="SimAss.CTL",append=TRUE)
 cat(InitSmooth,"\n",file="SimAss.CTL",append=TRUE) 
 cat("#Fishery Independent data available?","\n",file="SimAss.CTL",append=TRUE)
 cat(FisheryIndepenDat,"\n",file="SimAss.CTL",append=TRUE)
 cat("#Length at age standard deviation","\n",file="SimAss.CTL",append=TRUE)
 cat(GrowthSDn[length(steepnessN)],"\n",file="SimAss.CTL",append=TRUE)
 cat("#smoothness on fishing mortality penalty","\n",file="SimAss.CTL",append=TRUE)
 cat(FmortPen,"\n",file="SimAss.CTL",append=TRUE)
 cat("#smoothness on recruit penalty","\n",file="SimAss.CTL",append=TRUE)
 cat(RecruitPen,"\n",file="SimAss.CTL",append=TRUE)
 cat("#smoothness on time-varying natural mortality penalty","\n",file="SimAss.CTL",append=TRUE)
 cat(Mpenalty,"\n",file="SimAss.CTL",append=TRUE)
 cat("#smoothness on time-varying growth penalty","\n",file="SimAss.CTL",append=TRUE)
 cat(Growthpenalty,"\n",file="SimAss.CTL",append=TRUE)
 cat("#smoothness on time-varying selelctivity penalty","\n",file="SimAss.CTL",append=TRUE)
 cat(SelPenalty,"\n",file="SimAss.CTL",append=TRUE)

 cat("#harvest control rule selected","\n",file="SimAss.CTL",append=TRUE)
 cat(HarvestControlN,"\n",file="SimAss.CTL",append=TRUE)
 cat("#Constantcatch if harvest control ==2, set catch level","\n",file="SimAss.CTL",append=TRUE)
 cat(ConstantCatchN,"\n",file="SimAss.CTL",append=TRUE)
 cat("#Constantcatch if harvest control ==3, set fmort level","\n",file="SimAss.CTL",append=TRUE)
 cat(ConstantFn,"\n",file="SimAss.CTL",append=TRUE)
 cat("#40 10 parameters for harvest control 4","\n",file="SimAss.CTL",append=TRUE)
 cat(HCalphaN,"\n",file="SimAss.CTL",append=TRUE)
 cat("#40 10 parameters for harvest control 4","\n",file="SimAss.CTL",append=TRUE)
 cat(HCbetaN,"\n",file="SimAss.CTL",append=TRUE)

#==write .PIN file to get initial values close for estimation
 cat("# Simulated assessment pin file","\n",file="SimAss.PIN")
 cat("#","\n",file="SimAss.PIN",append=TRUE)

 InputSurvSel50N<-selAtAgeFunc(surv50n,VonKn,LinfN,t0n)
 InputSurvSel95N<-selAtAgeFunc(surv95n,VonKn,LinfN,t0n)

 cat("# srv_sel50","\n",file="SimAss.PIN",append=TRUE)
 cat(mean(InputSurvSel50N,na.rm=T)*rnorm(1,1,InitValSig),"\n",file="SimAss.PIN",append=TRUE)
 cat("# srv_sel95","\n",file="SimAss.PIN",append=TRUE)
 cat(mean(InputSurvSel95N,na.rm=T)*rnorm(1,1,InitValSig),"\n",file="SimAss.PIN",append=TRUE)
 cat("# stNatLen","\n",file="SimAss.PIN",append=TRUE)
 cat(log(VirInitN)*rnorm(length(log(VirInitN)),1,InitValSig),"\n",file="SimAss.PIN",append=TRUE)
 cat("# log_avg_fmort_dir","\n",file="SimAss.PIN",append=TRUE)
 cat(log(mean(trueFmortN,na.rm=T))*rnorm(1,1,InitValSig),"\n",file="SimAss.PIN",append=TRUE)
 cat("# fmort_dir_dev","\n",file="SimAss.PIN",append=TRUE)
 cat(rep(0,y-start_assessment),"\n",file="SimAss.PIN",append=TRUE)
 cat("# mean_log_rec","\n",file="SimAss.PIN",append=TRUE)
 cat(log(mean(trueRecN,na.rm=T))*rnorm(1,1,InitValSig),"\n",file="SimAss.PIN",append=TRUE)
 cat("# rec_dev","\n",file="SimAss.PIN",append=TRUE)
 cat(rep(0,y-start_assessment),"\n",file="SimAss.PIN",append=TRUE)
 cat("# log_avg_NatM","\n",file="SimAss.PIN",append=TRUE)
 if(EstM > 0) {
 cat(log(mean(NatMn,na.rm=T))*rnorm(1,1,InitValSig),"\n",file="SimAss.PIN",append=TRUE)
 }
 if(EstM <= 0) {
 cat(log(mean(NatMn,na.rm=T)),"\n",file="SimAss.PIN",append=TRUE) 
 }
 cat("# NatM_dev","\n",file="SimAss.PIN",append=TRUE)
 cat(rep(0,y-start_assessment),"\n",file="SimAss.PIN",append=TRUE)
 cat("# log_avg_GrowthK","\n",file="SimAss.PIN",append=TRUE)
 if(EstGrowthK > 0) {
 cat(log(mean(VonKn,na.rm=T))*rnorm(1,1,InitValSig),"\n",file="SimAss.PIN",append=TRUE)
 }
 if(EstGrowthK <= 0) {
 cat(log(mean(VonKn,na.rm=T)),"\n",file="SimAss.PIN",append=TRUE)
 }
 cat("# log_avg_Linf","\n",file="SimAss.PIN",append=TRUE)
 if(EstLinf > 0) {
 cat(log(mean(LinfN,na.rm=T))*rnorm(1,1,InitValSig),"\n",file="SimAss.PIN",append=TRUE)
 }
 if(EstLinf < 0) {
 cat(log(mean(LinfN,na.rm=T)),"\n",file="SimAss.PIN",append=TRUE)
 }
 cat("# GrowthK_dev","\n",file="SimAss.PIN",append=TRUE)
 cat(rep(0,y-start_assessment),"\n",file="SimAss.PIN",append=TRUE)
 cat("# Linf_dev","\n",file="SimAss.PIN",append=TRUE)
 cat(rep(0,y-start_assessment),"\n",file="SimAss.PIN",append=TRUE)

 InputSel50N<-selAtAgeFunc(sel50n,VonKn,LinfN,t0n)
 InputSel95N<-selAtAgeFunc(sel95n,VonKn,LinfN,t0n)

 cat("# SelPars50","\n",file="SimAss.PIN",append=TRUE)
 cat(mean(InputSel50N,na.rm=T)*rnorm(1,1,InitValSig),"\n",file="SimAss.PIN",append=TRUE)
 cat("# SelPars95","\n",file="SimAss.PIN",append=TRUE)
 cat(mean(InputSel95N,na.rm=T)*rnorm(1,1,InitValSig),"\n",file="SimAss.PIN",append=TRUE)

 cat("# SelPars_dev50","\n",file="SimAss.PIN",append=TRUE)
 cat(rep(0,y-start_assessment),"\n",file="SimAss.PIN",append=TRUE)
 cat("# SelPars_dev95","\n",file="SimAss.PIN",append=TRUE)
 cat(rep(0,y-start_assessment),"\n",file="SimAss.PIN",append=TRUE)
 	
 #=======================================
 #==.DAT file 
 #======================================
 file.create("SimAss.DAT")
 cat("#Simulated generalized assessment","\n",file="SimAss.DAT")
 cat("#","\n",file="SimAss.DAT",append=TRUE)
 cat("#start/end year","\n",file="SimAss.DAT",append=TRUE)
 cat(c(start_assessment,(y-1)),"\n",file="SimAss.DAT",append=TRUE)

 #==aggregate the length frequencies from the two areas
 cat("#number of ages","\n",file="SimAss.DAT",append=TRUE)
 cat(MaxAge,"\n",file="SimAss.DAT",append=TRUE)
 cat("#ages","\n",file="SimAss.DAT",append=TRUE)
 cat(seq(1,MaxAge),"\n",file="SimAss.DAT",append=TRUE)
 cat("#Number of length bins","\n",file="SimAss.DAT",append=TRUE)
 cat(LengthBinN,"\n",file="SimAss.DAT",append=TRUE)

 #==aggregate the index of abundance from two areas
 cat("#CPUE years","\n",file="SimAss.DAT",append=TRUE)
 cat((y-start_assessment),"\n",file="SimAss.DAT",append=TRUE)
 cat("#years for CPUE","\n",file="SimAss.DAT",append=TRUE)
 cat(seq(start_assessment,(y-1)),"\n",file="SimAss.DAT",append=TRUE)

 cat("#survey years","\n",file="SimAss.DAT",append=TRUE)
 cat((y-start_assessment),"\n",file="SimAss.DAT",append=TRUE)
 cat("#years for survey","\n",file="SimAss.DAT",append=TRUE)
 cat(seq(start_assessment,(y-1)),"\n",file="SimAss.DAT",append=TRUE)

 cat("#length freq from catch years","\n",file="SimAss.DAT",append=TRUE)
 cat((y-start_assessment),"\n",file="SimAss.DAT",append=TRUE)
 cat("#years for length freq catch","\n",file="SimAss.DAT",append=TRUE)
 cat(seq(start_assessment,(y-1)),"\n",file="SimAss.DAT",append=TRUE)
 cat("#catch length sample sizes","\n",file="SimAss.DAT",append=TRUE)
 cat(LenSampleN[seq(start_assessment,(y-1))],"\n",file="SimAss.DAT",append=TRUE) 

 cat("#length freq from survey years","\n",file="SimAss.DAT",append=TRUE)
 cat((y-start_assessment),"\n",file="SimAss.DAT",append=TRUE)
 cat("#years for length freq survey","\n",file="SimAss.DAT",append=TRUE)
 cat(seq(start_assessment,(y-1)),"\n",file="SimAss.DAT",append=TRUE)
 cat("#survey sample sizes","\n",file="SimAss.DAT",append=TRUE)
 cat(LenSampleN[seq(start_assessment,(y-1))],"\n",file="SimAss.DAT",append=TRUE) 

 cat("#total survey biomass","\n",file="SimAss.DAT",append=TRUE) 
 cat(SurvDataN[start_assessment:(y-1)],"\n",file="SimAss.DAT",append=TRUE)
 cat("#survey CV","\n",file="SimAss.DAT",append=TRUE)
 cat(IndexCVn[seq(start_assessment,(y-1))],"\n",file="SimAss.DAT",append=TRUE)

 cat("#total CPUE","\n",file="SimAss.DAT",append=TRUE) 
 cat(CPUEDataN[start_assessment:(y-1)],"\n",file="SimAss.DAT",append=TRUE)
 cat("#CPUE CV","\n",file="SimAss.DAT",append=TRUE)
 cat(IndexCVn[seq(start_assessment,(y-1))],"\n",file="SimAss.DAT",append=TRUE)

 cat("#total catch biomass","\n",file="SimAss.DAT",append=TRUE)
 cat(CatchDataN[start_assessment:(y-1)],"\n",file="SimAss.DAT",append=TRUE)
 cat("#catch CV","\n",file="SimAss.DAT",append=TRUE)
 cat(CatchCVn[seq(start_assessment,(y-1))],"\n",file="SimAss.DAT",append=TRUE)
 cat("#Length Bins","\n",file="SimAss.DAT",append=TRUE)
 cat(LengthBinsMid,"\n",file="SimAss.DAT",append=TRUE)

projCatLenFreqN[is.na(projCatLenFreqN)]<-0
projSurvLenFreqN[is.na(projSurvLenFreqN)]<-0

 cat("#catch length counts","\n",file="SimAss.DAT",append=TRUE)
 for(i in start_assessment:(y-1))
  cat(projCatLenFreqN[i,,z],"\n",file="SimAss.DAT",append=TRUE)

 cat("#survey length counts","\n",file="SimAss.DAT",append=TRUE)
 for(i in start_assessment:(y-1))
  cat(projSurvLenFreqN[i,,z],"\n",file="SimAss.DAT",append=TRUE)

 #==run the code
 if(os == "windows") {
  shell(SimAssComm)
 
 #==delete the .exe
  shell(paste("del", SimAssExec, sep = " "))
  }

  if(os != "windows") {
   system(SimAssComm)
   system(paste("rm", SimAssExec, sep = " "))
  }

#==================================
# PLOT THE ASSESSMENT OUTPUT (if instructed)
#==================================
REP<-readLines(file.path(MSEdir,IndSimFolder,"simass.rep"))
CTL<-readLines(file.path(MSEdir,IndSimFolder,"simass.CTL"))
DAT<-readLines(file.path(MSEdir,IndSimFolder,"simass.DAT"))
TRU<-readLines(file.path(MSEdir,IndSimFolder,"TrueQuantities.DAT"))
data2<-scan(file.path(MSEdir,IndSimFolder,"simass.PAR"),what="character")

 if(EstimationPlots==1 & z==Nsim & y == SimYear)
  AgeAssPlot(REP,CTL,DAT,TRU,data2,MSEdir,CreateFolderName)

temp<-grep("OFL",REP)[2]
OFL<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))

 #==calculate the TAC based on a SRR or a proxy, depending on data
 TAC[z,y]	<-OFL
 }


 #==================================================================================
 #==find the F that removes the TAC given the fishery selectivity===================
 #==this TAC is not necessarily the 'true' TAC--there is error in the biomass and FMSY
 #==================================================================================

 minExp<-0.00001
 maxExp<-0.99999
 inputFs	<-0  # NEED TO BE CHANGED TO THE INPUT F FOR POP 2 OR SOME OUTPUT OF HCR

 for(p in 1:30)
 {
 tempExp	<-(maxExp+minExp)/2
 tempF	<- -log(1-tempExp)
 tempFn	<-tempF
 tempFs	<-inputFs

   finCatLenFn	<-((tempFn*vulnN[y,])/(vulnN[y,]*tempFn+NatMn[y])) * (1-exp(-((tempFn*vulnN[y,])+NatMn[y]))) * projNn[y-1,,z]
   finCatFn		<-sum(finCatLenFn*WeightAtAgeN[y,])

   finCatLenFs	<-((tempFs*vulnS[y,])/(vulnS[y,]*tempFs+NatMs[y])) * (1-exp(-((tempFs*vulnS[y,])+NatMs[y]))) * projNs[y-1,,z]
   finCatFs		<-sum(finCatLenFs*WeightAtAgeS[y,])

 if(finCatFn<TAC[z,y])
  minExp	<-tempExp
 if(finCatFn>TAC[z,y])
  maxExp	<-tempExp
 initExp	<-tempExp
 }
 ApplyFn	<- -log(1-tempExp)
 ApplyFs	<- inputFs
 trueFmortN[z,y]	<-ApplyFn

 #==calculate true OFL
 #==find FOFL given spawning biomass
  tempSpBio			<-EggsN
  FutMort			<-trueF35[y]
  if(tempSpBio<trueB35[z,y-1])
  {
   FutMort<-0
   if(tempSpBio>HCbetaN*trueB35[z,y-1])
    FutMort = trueF35[y]*(tempSpBio/trueB35[z,y-1]-HCalphaN)/(1-HCalphaN)
   }
  #==find the catch for the FOFL
   trueCatchAtAge		<-((vulnN[y,]*FutMort)/(vulnN[y,]*FutMort+NatMn[y])) * (1-exp(-(vulnN[y,]*FutMort+NatMn[y]))) * projNn[y-1,,z]
   trueCatch		<-sum(trueCatchAtAge*WeightAtAgeN[y,])
   trueOFL[z,y]		<-trueCatch

 #===UPDATE DATA FOR ASSESSMENT====================================
   projCatchAtAgeN[y,,z]	<-((ApplyFn*vulnN[y,])/(vulnN[y,]*ApplyFn+NatMn[y])) * (1-exp(-((ApplyFn*vulnN[y,])+NatMn[y]))) * projNn[y-1,,z]
   projCatchN[z,y]		<-sum(projCatchAtAgeN[y,,z]*WeightAtAgeN[y,])

   projCatchAtAgeS[y,,z]	<-((ApplyFs*vulnS[y,])/(vulnS[y,]*ApplyFs+NatMs[y])) * (1-exp(-((ApplyFs*vulnS[y,])+NatMs[y]))) * projNs[y-1,,z]
   projCatchS[z,y]		<-sum(projCatchAtAgeS[y,,z]*WeightAtAgeS[y,])

   #==UPDATE CATCH DATA AND INDEX (survey and CPUE) DATA
   projSSBn[z,y]	<-sum(projNn[y-1,,z]*matureN[y,]*WeightAtAgeN[y,])
   projSSBs[z,y]	<-sum(projNs[y-1,,z]*matureS[y,]*WeightAtAgeS[y,])
   projExpBn[z,y]	<-sum(projNn[y-1,,z]*vulnN[y,]*WeightAtAgeN[y,])
   projExpBs[z,y]	<-sum(projNs[y-1,,z]*vulnS[y,]*WeightAtAgeS[y,])
   projSurvN[z,y]	<-sum(projNn[y-1,,z]*survSelN[y,]*WeightAtAgeN[y,])
   projSurvS[z,y]	<-sum(projNs[y-1,,z]*survSelS[y,]*WeightAtAgeS[y,])

   trueCatchN[z,y]	<-projCatchN[z,y]
   CatchAssessN[z,y]	<-projCatchN[z,y]*exp(CatchErrorN[z,y]-(CatchCVn[y]^2/2))
   CatchAssessS[z,y]	<-projCatchS[z,y]*exp(CatchErrorS[z,y]-(CatchCVs[y]^2/2))

   trueCPUEindN[z,y]	<-projExpBn[z,y]
   CPUEAssessN[z,y]	<-projExpBn[z,y]*exp(CPUEErrorN[z,y]-(IndexCVn[y]^2/2))
   CPUEAssessS[z,y]	<-projExpBs[z,y]*exp(CPUEErrorS[z,y]-(IndexCVs[y]^2/2))

   trueSurvIndN[z,y]	<-projSurvN[z,y]+projSurvS[z,y]
   SurvAssessN[z,y]	<-projSurvN[z,y]*exp(SurvErrorN[z,y]-(IndexCVn[y]^2/2))
   SurvAssessS[z,y]	<-projSurvS[z,y]*exp(SurvErrorS[z,y]-(IndexCVs[y]^2/2))
   
   trueSpbioN[z,y]	<-projSSBn[z,y]

  #==UPDATE CATCH LENGTH FREQS====
    for(w in 1:ncol(LenAtAgeN))
	 {
	  probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeN[y,w],sd=GrowthSDn[y])
	  ProbN<-probtemp/sum(probtemp)
	  tempCatAtLenN[w,,z]<-projCatchAtAgeN[y,w,z]*ProbN
	  probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeS[y,w],sd=GrowthSDs[y])
	  ProbS<-probtemp/sum(probtemp)
	  tempCatAtLenS[w,,z]<-projCatchAtAgeS[y,w,z]*ProbS
	 }
	 projCatLenFreqN[y,,z]<-apply(tempCatAtLenN[,,z],2,sum)
	 projCatLenFreqS[y,,z]<-apply(tempCatAtLenS[,,z],2,sum)

  #==============================================================
  # calculate the survey proportion at age (length) for assessment
  #==============================================================
	 #==make length frequencies for s==
   for(w in 1:ncol(LenAtAgeS))
    {
	probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeN[y,w],sd=GrowthSDn[y])
	ProbN<-probtemp/sum(probtemp)
	tempSurvAtLenN[w,,z]<-projNn[y-1,w,z]*survSelN[y,w]*ProbN
	probtemp<-dnorm(LengthBinsMid,mean=LenAtAgeS[y,w],sd=GrowthSDs[y])
	ProbS<-probtemp/sum(probtemp)
	tempSurvAtLenS[w,,z]<-projNs[y-1,w,z]*survSelS[y,w]*ProbS
     }

    projSurvLenFreqN[y,,z]<-apply(tempSurvAtLenN[,,z],2,sum)
    projSurvLenFreqS[y,,z]<-apply(tempSurvAtLenS[,,z],2,sum)

   # sum(projSurvLenFreqN[y,,z])
   # sum(projNn[y,,z])

 #==apply F that would take a calculated TAC====================================================
 #==update population dynamics
  for (i in 2:(MaxAge-1))
  {
   projNn[y,i,z]		<-projNn[y-1,i-1,z]*exp(-ApplyFn*vulnN[y,i-1])*exp(-NatMn[y])
   projNs[y,i,z]		<-projNs[y-1,i-1,z]*exp(-ApplyFs*vulnS[y,i-1])*exp(-NatMs[y])
  }
   projNn[y,MaxAge,z]	<-(projNn[y-1,(MaxAge-1),z])*exp(-ApplyFn*vulnN[y,MaxAge])*exp(-NatMn[y])+ projNn[y-1,MaxAge,z]*exp(-ApplyFn*vulnN[y,MaxAge])*exp(-NatMn[y])
   projNs[y,MaxAge,z]	<-(projNs[y-1,(MaxAge-1),z])*exp(-ApplyFs*vulnS[y,MaxAge])*exp(-NatMs[y])+ projNs[y-1,MaxAge,z]*exp(-ApplyFn*vulnS[y,MaxAge])*exp(-NatMs[y])

   moveFromN		<-projNn[y,,z]*MovementN[y,]
   moveFromS		<-projNs[y,,z]*MovementS[y,]

   projNn[y,,z]		<-projNn[y,,z]-moveFromN+moveFromS
   projNs[y,,z]		<-projNs[y,,z]-moveFromS+moveFromN

   EggsN			<-sum(projNn[y-1,,z]*matureN[y,]*WeightAtAgeN[y,])
   EggsS			<-sum(projNs[y-1,,z]*matureS[y,]*WeightAtAgeS[y,])

   projNn[y,1,z]		<-Recruitment(EggsIN=EggsN,steepnessIN=steepnessN[y],RzeroIN=RzeroN[y],RecErrIN=RecErrN[z,y],recType="BH",NatMin=NatMn[y],
							vulnIN=vulnN[y,],matureIN=matureN[y,],weightIN=WeightAtAgeN[y,],LenAtAgeIN=LenAtAgeN[y,],MaxAge=MaxAge,sigmaRin=sigmaRn[y])
   projNs[y,1,z]		<-Recruitment(EggsIN=EggsS,steepnessIN=steepnessS[y],RzeroIN=RzeroS[y],RecErrIN=RecErrS[z,y],recType="BH",NatMin=NatMs[y],
							vulnIN=vulnS[y,],matureIN=matureS[y,],weightIN=WeightAtAgeS[y,],LenAtAgeIN=LenAtAgeS[y,],MaxAge=MaxAge,sigmaRin=sigmaRs[y])
   trueRecN[z,y]		<-projNn[y,1,z]

  if(!silent) cat(paste("Year ",y," of ",SimYear," in simulation ",z," of ",Nsim,"\n",sep=""))
  
  trueB35[z,y]  <-trueSBPR35[y]*mean(trueRecN[z,start_assessment:SimYear],na.rm=T)
 } # end y
} # end z

} # end of GeMS


