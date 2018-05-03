#' More plots.. currently commented out from \code{\link{AgeStructureComp}}
#'
#' @param CTLName Name of CTL file
#' @param out List output from \code{\link{ReadCTLfile}}
#' @param MSEdir Directory containing CTL files
#'
#' @return The path to an GeMS Age-Structured Model binary.
#'
#' @export
#' @examples
#' \dontrun{
#' GeMSplots("Cod_Base_CTL",ReadCTLfile("Cod_Base_CTL.csv"),"~/GeMS/Examples/Cod_1_Production")
#' }
#' 

GeMSplots<-function(CTLName,out,MSEdir)
{
Nsim			<-out$OM$Nsim			# number of simulations to do in the MSE
SimYear		<-out$OM$SimYear			# total number of years in simulation
InitYear		<-out$OM$InitYear			# year in which MSE starts (i.e. the number of years of data available)
AssessmentType	<-out$OM$AssessmentType		# this can be "Production" or "Non-spatial age" ("Spatial age" and "Time-varying production" in the works)
MaxAge		<-out$OM$MaxAge

if(AssessmentType==1)
{
 PROD<-readLines(file.path(MSEdir,CTLName,"ProdOuts.csv"))

temp<-grep("true Catch",PROD)[1]
TrueCatchOut<-matrix(as.numeric(unlist(strsplit(PROD[(temp+1):(temp+Nsim)],split=" "))),nrow=Nsim,byrow=T)
temp<-grep("true fishing mortality",PROD)[1]
TrueFmortOut<-matrix(as.numeric(unlist(strsplit(PROD[(temp+1):(temp+Nsim)],split=" "))),nrow=Nsim,byrow=T)
temp<-grep("true spawning biomass",PROD)[1]
TrueSpBioOut<-matrix(as.numeric(unlist(strsplit(PROD[(temp+1):(temp+Nsim)],split=" "))),nrow=Nsim,byrow=T)
temp<-grep("true CPUE",PROD)[1]
TrueCPUEOut<-matrix(as.numeric(unlist(strsplit(PROD[(temp+1):(temp+Nsim)],split=" "))),nrow=Nsim,byrow=T)
temp<-grep("true total allowable catch",PROD)[1]
TrueTACOut<-matrix(as.numeric(unlist(strsplit(PROD[(temp+1):(temp+Nsim)],split=" "))),nrow=Nsim,byrow=T)

temp<-grep("true BMSY",PROD)[1]
TrueBMSYOut<-matrix(as.numeric(unlist(strsplit(PROD[(temp+1):(temp+Nsim)],split=" "))),nrow=Nsim,byrow=T)
temp<-grep("true FMSY",PROD)[1]
TrueFMSYOut<-matrix(as.numeric(unlist(strsplit(PROD[(temp+1):(temp+Nsim)],split=" "))),nrow=Nsim,byrow=T)

temp<-grep("est cpue ind",PROD)[1]
EstCPUEOut<-matrix(as.numeric(unlist(strsplit(PROD[(temp+1):(temp+Nsim)],split=" "))),nrow=Nsim,byrow=T)
temp<-grep("est spawning biomass by year",PROD)[1]
EstCurbioOut<-matrix(as.numeric(unlist(strsplit(PROD[(temp+1):(temp+Nsim)],split=" "))),nrow=Nsim,byrow=T)
temp<-grep("est total allowable catch",PROD)[1]
EstTACOut<-matrix(as.numeric(unlist(strsplit(PROD[(temp+1):(temp+Nsim)],split=" "))),nrow=Nsim,byrow=T)
temp<-grep("est BMSY",PROD)[1]
EstBMSYOut<-matrix(as.numeric(unlist(strsplit(PROD[(temp+1):(temp+Nsim)],split=" "))),nrow=Nsim,byrow=T)
temp<-grep("est FMSY",PROD)[1]
EstFMSYOut<-matrix(as.numeric(unlist(strsplit(PROD[(temp+1):(temp+Nsim)],split=" "))),nrow=Nsim,byrow=T)

temp<-grep("data cpue",PROD)[1]
DateCPUEout<-matrix(as.numeric(unlist(strsplit(PROD[(temp+1):(temp+Nsim)],split=" "))),nrow=Nsim,byrow=T)

pdf(file.path(MSEdir,CTLName,"RelativeErrors.pdf"),height=3,width=4)
PolygonPlots(Truth=TrueCPUEOut,Estimated=EstCPUEOut,Observed=DateCPUEout)
dev.off()

#==plots
REtac<-(EstTACOut-TrueTACOut)/TrueTACOut
REbmsy<-(EstBMSYOut-TrueBMSYOut[1])/TrueBMSYOut[1]
REfmsy<-(EstFMSYOut-TrueFMSYOut[1])/TrueFMSYOut[1]
REbio<-(EstCurbioOut-TrueSpBioOut)/TrueSpBioOut
LostYield<-TrueTACOut-EstTACOut
PercLostYield<-LostYield/TrueTACOut
PercOverfishedTrue<-TrueSpBioOut/TrueBMSYOut[1]
PercOverfishedEst<-EstCurbioOut/TrueBMSYOut[1]

pdf(file.path(MSEdir,CTLName,"RelativeErrors.pdf"),height=8,width=4)
par(mfrow=c(4,1),mar=c(0.1,4,.1,1),oma=c(3,0,3,0))
MyBoxPlot(REtac[,InitYear:SimYear],"Total allowable catch")
MyBoxPlot(REbmsy[,InitYear:SimYear],"BMSY")
MyBoxPlot(REfmsy[,InitYear:SimYear],"FMSY")
MyBoxPlot(REbio[,InitYear:SimYear],"Estimated biomass",bottomplot=1)
mtext(outer=T,side=3,"Relative error")

pdf(file.path(MSEdir,CTLName,"YieldStats.pdf"),height=8,width=4)
par(mfrow=c(4,1),mar=c(0.1,4,.1,1),oma=c(3,0,0,0))
MyBoxPlot(LostYield[,InitYear:SimYear],"Lost yield")
MyBoxPlot(PercLostYield[,InitYear:SimYear],"Lost yield/yield")
MyBoxPlot(PercOverfishedTrue[,InitYear:SimYear],"True B/BMSY",refVal=1)
MyBoxPlot(PercOverfishedEst[,InitYear:SimYear],"Estimated B/BMSY",bottomplot=1,refVal=1)

}

if(AssessmentType==2)
{


#==predicted derived quantities storage
predSpBio		<-matrix(ncol=SimYear-1,nrow=Nsim)
predSurvBio		<-matrix(ncol=SimYear-1,nrow=Nsim)
obsSurvBio		<-matrix(ncol=SimYear-1,nrow=Nsim)
predCpueBio 	<-matrix(ncol=SimYear-1,nrow=Nsim)
obsCpueIndex 	<-matrix(ncol=SimYear-1,nrow=Nsim)
predCatchBio	<-matrix(ncol=SimYear-1,nrow=Nsim)
obsCatchBio		<-matrix(ncol=SimYear-1,nrow=Nsim)
obsSurlen		<-array(dim=c(SimYear-1,MaxAge,Nsim))
predSurlen		<-array(dim=c(SimYear-1,MaxAge,Nsim))
obsCatchlen		<-array(dim=c(SimYear-1,MaxAge,Nsim))
predCatchlen	<-array(dim=c(SimYear-1,MaxAge,Nsim))
trueRecPlot		<-matrix(ncol=SimYear-1,nrow=Nsim)
trueCatchPlot	<-matrix(ncol=SimYear-1,nrow=Nsim)
trueSpbioPlot	<-matrix(ncol=SimYear-1,nrow=Nsim)
trueCPUEindPlot	<-matrix(ncol=SimYear-1,nrow=Nsim)
trueSurvIndPlot	<-matrix(ncol=SimYear-1,nrow=Nsim)
trueFmortPlot	<-matrix(ncol=SimYear-1,nrow=Nsim)

trueB35Plot		<-rep(0,Nsim)
trueBMSYPlot	<-rep(0,Nsim)
trueF35Plot		<-rep(0,Nsim)
trueFMSYPlot	<-rep(0,Nsim)

predB35Plot		<-matrix(ncol=SimYear-1,nrow=Nsim)
predBMSYPlot	<-matrix(ncol=SimYear-1,nrow=Nsim)
predF35Plot		<-matrix(ncol=SimYear-1,nrow=Nsim)
predFMSYPlot	<-matrix(ncol=SimYear-1,nrow=Nsim)

predB35PlotEnd	<-rep(0,Nsim)
predBMSYPlotEnd	<-rep(0,Nsim)
predF35PlotEnd	<-rep(0,Nsim)
predFMSYPlotEnd	<-rep(0,Nsim)


#===parameters===
meanREC		<-rep(0,Nsim)
recDevs		<-matrix(ncol=SimYear-1,nrow=Nsim)
TotRec		<-matrix(ncol=SimYear-1,nrow=Nsim)
logFdir		<-rep(0,Nsim)
fmort_dir_dev	<-matrix(ncol=SimYear-1,nrow=Nsim)
TotFdir		<-matrix(ncol=SimYear-1,nrow=Nsim)
sel50surv		<-rep(0,Nsim)
sel95surv		<-rep(0,Nsim)
sel50fish		<-rep(0,Nsim)
sel95fish		<-rep(0,Nsim)

OutFolder<-file.path(MSEdir,CTLName)
pdf(file.path(OutFolder,"AgeStrucurePlots.pdf"))
#==cycle through final .rep files, pull out data
for(n in 1:Nsim)
{
IndSimFolder<-file.path(MSEdir,CTLName,n,SimYear)

REP<-readLines(file.path(IndSimFolder,"simass.rep"))
CTL<-readLines(file.path(IndSimFolder,"simass.CTL"))
DAT<-readLines(file.path(IndSimFolder,"simass.DAT"))
TRU<-readLines(file.path(IndSimFolder,"TrueQuantities.DAT"))
data2<-scan(file.path(IndSimFolder,"simass.PAR"),what="character")

if(n==1)
{
temp<-grep("survey years",DAT)[1]
yearsDat<-as.numeric(unlist(strsplit(DAT[temp+1],split=" ")))
temp<-grep("number of ages",DAT)[1]
maxAge<-as.numeric(unlist(strsplit(DAT[temp+1],split=" ")))
AgeBin<-seq(1,maxAge)
temp<-grep("Length Bin",DAT)[1]
LengthBin<-as.numeric(unlist(strsplit(DAT[temp+1],split=" ")))
}

temp<-grep("recruitment",TRU)
trueRecPlot[n,]		<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))[1:(yearsDat)]
temp<-grep("Catch",TRU)
trueCatchPlot[n,]		<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))[1:(yearsDat)]
temp<-grep("fishing mortality",TRU)
trueFmortPlot[n,]	<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))[1:(yearsDat)]
temp<-grep("spawning biomass",TRU)
trueSpbioPlot[n,]		<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))[1:(yearsDat)]
temp<-grep("cpue index",TRU)
trueCPUEindPlot[n,]	<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))[1:(yearsDat)]
temp<-grep("survey index",TRU)
trueSurvIndPlot[n,]	<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))[1:(yearsDat)]

#temp<-grep("B35",TRU)
#trueB35Plot[n]	<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))
#temp<-grep("F35",TRU)
#trueF35Plot[n]	<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))
temp<-grep("BMSY",TRU)
trueBMSYPlot[n]	<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))
temp<-grep("FMSY",TRU)
trueFMSYPlot[n]	<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))

temp<-grep("spawning biomass",REP)
predSpBio[n,]	<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+1)]
temp<-grep("B35",REP)
predB35PlotEnd[n]	<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))
temp<-grep("F35",REP)[2]
predF35PlotEnd[n]	<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))

temp<-grep("pred survey bio",REP)
predSurvBio[n,]	<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+1)]
temp<-grep("obs survey bio",REP)
obsSurvBio[n,]	<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+1)]

temp<-grep("pred cpue bio",REP)
predCpueBio[n,]	<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+1)]
temp<-grep("obs cpue bio",REP)
obsCpueIndex[n,] 	<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+1)]

temp<-grep("pred catch bio",REP)
predCatchBio[n,]	<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+1)]
temp<-grep("obs catch bio",REP)
obsCatchBio[n,]	<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+1)]

temp<-grep("obs surv len freq",REP)
obsSurlen[,,n]	<-matrix(as.numeric(unlist(strsplit(REP[(temp+1):(temp+yearsDat)],split=" "))),nrow=yearsDat,byrow=T)[,2:(maxAge+1)]
temp<-grep("pred surv len freq",REP)
predSurlen[,,n]	<-matrix(as.numeric(unlist(strsplit(REP[(temp+1):(temp+yearsDat)],split=" "))),nrow=yearsDat,byrow=T)[,2:(maxAge+1)]

temp<-grep("obs catch len freq",REP)
obsCatchlen[,,n]	<-matrix(as.numeric(unlist(strsplit(REP[(temp+1):(temp+yearsDat)],split=" "))),nrow=yearsDat,byrow=T)[,2:(maxAge+1)]
temp<-grep("pred catch len freq",REP)
predCatchlen[,,n]	<-matrix(as.numeric(unlist(strsplit(REP[(temp+1):(temp+yearsDat)],split=" "))),nrow=yearsDat,byrow=T)[,2:(maxAge+1)]

#==================================================
#   Par file plots
#================================================
cnt<-grep("mean_log_rec",data2)
meanREC[n]<-as.numeric(data2[cnt+1])
cnt<-grep("rec_dev",data2)
recDevs[n,]<-as.numeric(data2[(cnt+1):(cnt+yearsDat)])
TotRec[n,]<-exp(meanREC[n]+recDevs[n,])

fmortYrs<-seq(1,yearsDat)
cnt<-grep("log_avg_fmort_dir",data2)
logFdir[n]<-as.numeric(data2[cnt+1])
cnt<-grep("fmort_dir_dev",data2)
fmort_dir_dev[n,]<-as.numeric(data2[(cnt+1):(cnt+length(fmortYrs))])
TotFdir[n,]<-exp(logFdir[n]+fmort_dir_dev[n,])
}

#==plot the fitted data
par(mfrow=c(3,1),mar=c(.1,3,.1,.1),oma=c(2,3,0,0))
PolygonPlots(Truth=trueCatchPlot,Estimated=predCatchBio,Observed=obsCatchBio,AddLegend=T,quantity="Catch")
PolygonPlots(Truth=trueCPUEindPlot,Estimated=predCpueBio,Observed=obsCpueIndex,quantity="CPUE")
PolygonPlots(Truth=trueSurvIndPlot,Estimated=predSurvBio,Observed=obsSurvBio,quantity="Survey")

#==plot the estimated processes
par(mfrow=c(2,1),mar=c(.1,3,.1,.1),oma=c(2,3,0,0))
PolygonPlots(Truth=trueRecPlot,Estimated=TotRec,AddLegend=T,quantity="Recruitment")
PolygonPlots(Truth=trueFmortPlot,Estimated=TotFdir,quantity="Fishing mortality")

#==SELECTIVITY==
# selAtAgeFunc(sel50n[1],VonKn[1],LinfN[1],t0n[1])
# selAtAgeFunc(sel95n[1],VonKn[1],LinfN[1],t0n[1])
# selAtAgeFunc(surv50n[1],VonKn[1],LinfN[1],t0n[1])
# selAtAgeFunc(surv95n[1],VonKn[1],LinfN[1],t0n[1])
# pull pars from OM
# translate to age with function
# calculate
# plot

#==plot the reference points
PolygonPlots(Truth=trueSpbioPlot,Estimated=predSpBio)
  plot(trueRecPlot~trueSpbioPlot,ylim=c(0,max(trueRecPlot,na.rm=T)))

#==errors in reference points
REbmsy<-(predB35PlotEnd-trueBMSYPlot)/trueBMSYPlot
REfmsy<-(predF35PlotEnd-trueFMSYPlot)/trueFMSYPlot
par(mfrow=c(2,1),mar=c(.1,3,.1,.1),oma=c(2,3,0,0))
boxplot(REbmsy)
boxplot(REfmsy)
## relative errors in everything
## summary table to compare between other scenarios
PlotSel<-0
if(PlotSel==1)
{
#dev.new()
par(mfrow=c(4,1),mar=c(1,4,1,1),oma=c(3,1,1,1))
#==plot estimated selectivities

cnt<-grep("srv_sel50",data2)
sel50surv<-as.numeric(data2[cnt+1])
cnt<-grep("srv_sel95",data2)
sel95surv<-as.numeric(data2[cnt+1])

cnt<-grep("SelPars50",data2)
sel50fish<-as.numeric(data2[cnt+1])
cnt<-grep("SelPars95",data2)
sel95fish<-as.numeric(data2[cnt+1])

trueSel50<-out$OM$sel50s
trueSel95<-out$OM$sel95s
trueSurv50<-out$OM$surv50s
trueSurv95<-out$OM$surv95s

SurvSel<-1/( 1 + exp( -1*log(19)*(AgeBin-sel50surv)/(sel95surv-sel50surv)))
FishSel<-1/( 1 + exp( -1*log(19)*(AgeBin-sel50fish)/(sel95fish-sel50fish)))

trueSurvSel<-1/( 1 + exp( -1*log(19)*(LengthBin-trueSurv50)/(trueSurv95-trueSurv50)))
trueFishSel<-1/( 1 + exp( -1*log(19)*(LengthBin-trueSel50)/(trueSel95-trueSel50)))

plot(trueSurvSel~LengthBin,type="l",xlab="Length",ylab="Selectivity")

lines(FishSel~AgeBin,lty=2,col=2)
lines(SurvSel~AgeBin,lty=2,col=1)
lines(trueFishSel~LengthBin,lty=1,col=2)
legend("bottomright",bty='n',col=c(1,1,2,2),lty=c(1,2,1,2),legend=c("True Survey","Est Survey","True Fishery","Est Fishery"))
}

#dev.new()
putIn<-ceiling(sqrt(yearsDat))
par(mfrow=c(putIn,putIn),mar=c(0.1,.1,.1,.1))
for(i in 2:yearsDat)
{
plot(obsSurlen[i,,1],axes=F)
lines(predSurlen[i,,1])
}
plot.new()
legend("center","Survey")
#dev.new()
putIn<-ceiling(sqrt(yearsDat))
par(mfrow=c(putIn,putIn),mar=c(0.1,.1,.1,.1))
for(i in 2:yearsDat)
{
plot(obsCatchlen[i,,1],axes=F)
lines(predCatchlen[i,,1])
}
plot.new()
legend("center","Catch")
dev.off()
}

}