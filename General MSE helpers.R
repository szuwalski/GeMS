###############################################################
SBPRfunc<-function(F,Rzero,NatM,vuln,mature,weight,LenAtAge)
{
#initial conditions
N<-matrix(ncol=MaxAge,nrow=50)
N[1,]<-initialN(Rzero,NatM,MaxAge)

for (j in 2:50)
{
 for (i in 2:(MaxAge-1))
  N[j,i]<-N[j-1,i-1]*exp(-F*vuln[i-1])*exp(-NatM)

 N[j,MaxAge]<-(N[j-1,MaxAge-1]+N[j-1,MaxAge])*exp(-F*vuln[MaxAge])*exp(-NatM)
 N[j,1]<-Rzero
}
SBP<-sum(N[50,]*mature*weight)
SBPR<-sum(N[50,]*mature*weight)/Rzero
YPR<-sum(N[50,]*(1-exp(-F*vuln)) *weight)
list(SBP,SBPR,YPR)
#return(SBPR)
}

###############################################################

Recruitment<-function(EggsIN,steepnessIN,RzeroIN,RecErrIN,recType,NatMin,vulnIN,matureIN,weightIN,LenAtAgeIN)
{
if(recType=="BH")
{
 Szero<-unlist(SBPRfunc(0,RzeroIN,NatMin,vulnIN,matureIN,weightIN,LenAtAgeIN)[1])	#since there are only a couple possible values, code could be sped up by changing this
 alpha<-(Szero/RzeroIN)*(1-steepnessIN)/(4*steepnessIN)
 beta<-(5*steepnessIN-1)/(4*steepnessIN*RzeroIN)
 Recruit<-(EggsIN/(alpha+beta*EggsIN))*exp(RecErrIN)
}
if(recType=="Rick")
{
 Szero<-unlist(SBPRfunc(0,RzeroIN,NatMin,vulnIN,matureIN,weightIN)[1])
 beta<-log(5*steepnessIN)/(0.8*Szero)
 alpha<-log(RzeroIN/Szero)+(beta*Szero)
 Recruit<-(EggsIN*exp(alpha-beta*EggsIN))*exp(RecErrIN)
}
return(Recruit)
}

###############################################################

initialN<-function(Rzero,NatM,inAge)
{
test		<-rep(0,100)
test[1]	<-Rzero[1]
for (i in 2:100)
 {test[i]	<-test[i-1]*exp(-NatM)}
virInit2	<-c(test[1:(inAge-1)],sum(test[inAge:100]))
return(virInit2)
}

###############################################################

BevHolt<-function(x)
{
alpha		<-abs(x[1])
beta		<-abs(x[2])
sigma		<-abs(x[3])
Pred		<-alpha*spbio/(1+(spbio/beta))
loglike	<-log(1/(sqrt(2*3.141598*sigma^2))) - ((log(Pred)-log(rec))^2)/(2*sigma^2) 
return(sum(-1*loglike))
}

###############################################################

Ricker<-function(x)
{
alpha		<-x[1]
beta		<-x[2]
sigma		<-abs(x[3])
Pred		<-spbio*exp(alpha-beta*spbio)
loglike	<-log(1/(sqrt(2*3.141598*sigma^2))) - ((log(Pred)-log(rec))^2)/(2*sigma^2) 
return(sum(-1*loglike))
}

###############################################################

StrLine<-function(x)
{
 Pred		<-x[1]
 sigma	<-x[2]
 loglike	<- log(1/(sqrt(2*3.141598*sigma^2))) - ((log(Pred)-log(rec))^2)/(2*sigma^2) 
 return(sum(-1*loglike))
}

###############################################################

StockRecCheck<-function(rec,spbio,plotSR=F)
{
Stx		<-srStarts(rec~spbio,type="BevertonHolt",param=2)
x		<-c(Stx$a,Stx$Rp,1)
outs		<-optim(x,BevHolt)
counter	<-1
while(outs$convergence!=0 & counter<20)
{
 x		<-outs$par
 outs		<-optim(x,BevHolt)
 counter	<-counter+1
}

Stx1<-srStarts(rec~spbio,type="Ricker",param=2)
x1		<-c(Stx1$a,Stx1$b,1)
outs1		<-optim(x1,Ricker)
x1		<-outs1$par
outs1		<-optim(x1,Ricker)
x1		<-outs1$par
outs1		<-optim(x1,Ricker)

x2		<-c(mean(rec),sd(rec))
outs2		<-optim(x2,StrLine)

AvgAIC	<-outs2$value
RickerAIC	<-outs1$value
BHaic		<-outs$value

AvgPar	<-outs2$par
RickerPars	<-outs1$par
BHpars	<-outs$par

#===check parameters
AICR		<-6+2*RickerAIC
AICB		<-6+2*BHaic
AICA		<-4+2*AvgAIC

#==report decision about stock recruit relationship
SRchoice	<-"Average"
returnPar	<-AvgPar

if(AICR<AICA-2 & AICR<AICB)
{
SRchoice	<-"Rick"
returnPar	<-RickerPars
}

if(AICB<AICA-2 & AICB<AICR)
{
SRchoice	<-"BevH"
returnPar	<-BHpars
}

if(plotSR==T)
 {
  predSpbio<-seq(1,max(spbio),max(spbio)/40)
  alpha<-abs(BHpars[1])
  beta<-abs(BHpars[2])
  plotPred<-alpha*predSpbio/(1+(predSpbio/beta))
  plot(rec~spbio,ylim=c(0,max(rec)),xlim=c(0,max(spbio)))
  lines(plotPred~predSpbio)

  alpha<-RickerPars[1]
  beta<-RickerPars[2]
  plotPred<-predSpbio*exp(alpha-beta*predSpbio)
  lines(plotPred~predSpbio,lty=2,col=2)


  abline(h=AvgPar[1],lty=4,col=4)

  legend("topright",bty='n',col=c(1,2,4),lty=c(1,2,4),legend=c(paste("BH ",round(AICB,2)),
   paste("Rick ", round(AICR,2)),paste("Avg ",round(AICA,2))))
 }
list(SRchoice,returnPar)
}

###############################################################

calcBHrec<-function(SRouts,spbioIn)
{
 alpha	<-abs(SRouts[[2]][1])
 beta		<-abs(SRouts[[2]][2])
 Pred		<-alpha*spbioIn/(1+(alpha*spbioIn/beta)) 
 return(Pred)
}

###############################################################

calcRICKrec<-function(SRouts,spbioIn)
{
 alpha	<-SRouts[[2]][1]
 beta		<-SRouts[[2]][2]
 Pred		<-spbioIn*exp(alpha-beta*spbioIn)
 return(Pred)
}

###############################################################

calcREF<-function(SRouts)
{
if(SRouts[[1]][1]!="Average")
{
 testExp	<-seq(0.001,1,0.01)
 recYield	<-rep(0,length(testExp))
 projYrs	<-100
 tempN	<-matrix(nrow=projYrs,ncol=maxAge)
 tempN[1,]	<-virInit
 tempYield	<-rep(0,length(testExp))
for(f in 1:length(testF))
{
 tempExp<-testExp[f]
 for(j in 2:projYrs)
 {
  for (i in 2:(maxAge-1))
   tempN[j,i]		<-tempN[j-1,i-1]*(1-tempExp*vuln[i-1])*survival
   tempN[j,maxAge]	<-(tempN[j-1,(maxAge-1)])*(1-tempExp*vuln[maxAge])*survival + tempN[j-1,maxAge]*(1-tempExp*vuln[maxAge])*survival
   spbioIn			<-sum(tempN[j-1,]*mature*weight)
   if(SRouts[[1]][1]=="BevH")
    tempN[j,1]		<-calcBHrec(SRouts,spbioIn) 
   if(SRouts[[1]][1]=="Rick")
    tempN[j,1]		<-calcRICKrec(SRouts,spbioIn) 
  }
   tempYield[f]		<-sum(tempN[j,]*tempExp*vuln*weight)
 }
}
calcRICKrec(SRouts,virBio)
if(SRouts[[1]][1]=="Average")
{
#STUFF
}

}

###############################################################
 
calcRemoval<-function( calcREFouts ) 
{
#==this should have options for how the control rule will change
 }


##############################################################
ProdMod<-function(x,CatchData,IndexData)
{
K<-abs(x[1])
r<-abs(x[2])
predBio<-rep(0,length(IndexData))
predBio[1]<-K
for(i in 2:length(CatchData))
{
 predBio[i]<-predBio[i-1]+r*predBio[i-1]*(1-predBio[i-1]/K)-CatchData[i]
}
SSQ<-sum((predBio-IndexData)^2,na.rm=T)
return(SSQ)
}

#################################################################
ProdModPlot<-function(x,CatchData,IndexData,plots=0)
{
K<-abs(x[1])
r<-abs(x[2])
predBio<-rep(0,length(IndexData)+1)
predBio[1]<-K
CatchData<-c(CatchData,0)
for(i in 2:length(CatchData))
{
 predBio[i]<-predBio[i-1]+r*predBio[i-1]*(1-predBio[i-1]/K)-CatchData[i]
}
if(plots==1)
{
plot(IndexData[1:(length(IndexData)-1)],ylim=c(0,max(IndexData,na.rm=T)),las=1,ylab="")
lines(predBio[1:(length(predBio)-1)])
lines(CatchData[1:(length(CatchData)-1)],lty=2,col=2)
}
return(predBio)
}

#################################################################
AgeAssPlot<-function(REP,CTL,DAT,TRU,data2)
{
temp<-grep("survey years",DAT)[1]
yearsDat<-as.numeric(unlist(strsplit(DAT[temp+1],split=" ")))
temp<-grep("number of ages",DAT)[1]
maxAge<-as.numeric(unlist(strsplit(DAT[temp+1],split=" ")))
AgeBin<-seq(1,maxAge)
temp<-grep("Length Bin",DAT)[1]
LengthBin<-as.numeric(unlist(strsplit(DAT[temp+1],split=" ")))

temp<-grep("recruitment",TRU)
trueRec<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))
temp<-grep("Catch",TRU)
trueCatch<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))[1:(yearsDat)]
temp<-grep("survey selectivity",TRU)
trueSurvSelPar<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))
temp<-grep("fishery selectivity",TRU)
trueFishSelPar<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))
temp<-grep("fishing mortality",TRU)
trueFmort<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))
temp<-grep("spawning biomass",TRU)
trueSpbio<-as.numeric(unlist(strsplit(TRU[temp+1],split=" ")))

temp<-grep("spawning biomass",REP)
predSpBio<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+1)]
temp<-grep("B35",REP)
B35<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))


temp<-grep("pred survey bio",REP)
predSurvBio<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+1)]
temp<-grep("obs survey bio",REP)
obsSurvBio<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+1)]

temp<-grep("pred cpue bio",REP)
predCpueBio <-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+1)]
temp<-grep("obs cpue bio",REP)
cpueIndex <-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+1)]

temp<-grep("pred catch bio",REP)
predCatchBio<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+1)]
temp<-grep("obs catch bio",REP)
obsCatchBio<-as.numeric(unlist(strsplit(REP[temp+1],split=" ")))[2:(yearsDat+1)]

temp<-grep("obs surv len freq",REP)
obsSurlen<-matrix(as.numeric(unlist(strsplit(REP[(temp+1):(temp+yearsDat)],split=" "))),nrow=yearsDat,byrow=T)[,2:(maxAge+1)]
temp<-grep("pred surv len freq",REP)
predSurlen<-matrix(as.numeric(unlist(strsplit(REP[(temp+1):(temp+yearsDat)],split=" "))),nrow=yearsDat,byrow=T)[,2:(maxAge+1)]

temp<-grep("obs catch len freq",REP)
obsCatchlen<-matrix(as.numeric(unlist(strsplit(REP[(temp+1):(temp+yearsDat)],split=" "))),nrow=yearsDat,byrow=T)[,2:(maxAge+1)]
temp<-grep("pred catch len freq",REP)
predCatchlen<-matrix(as.numeric(unlist(strsplit(REP[(temp+1):(temp+yearsDat)],split=" "))),nrow=yearsDat,byrow=T)[,2:(maxAge+1)]

#==================================================
#   Par file plots
#================================================
cnt<-grep("mean_log_rec",data2)
meanREC<-as.numeric(data2[cnt+1])
cnt<-grep("rec_dev",data2)
recDevs<-as.numeric(data2[(cnt+1):(cnt+yearsDat)])
TotRec<-exp(meanREC+recDevs)

fmortYrs<-seq(1,yearsDat)
cnt<-grep("log_avg_fmort_dir",data2)
logFdir<-as.numeric(data2[cnt+1])
cnt<-grep("fmort_dir_dev",data2)
fmort_dir_dev<-as.numeric(data2[(cnt+1):(cnt+length(fmortYrs))])
TotFdir<-exp(logFdir+fmort_dir_dev)

cnt<-grep("stNatLen",data2)
stNatLen<-as.numeric(data2[(cnt+1):(cnt+maxAge)])

dev.new()
par(mfrow=c(4,1),mar=c(1,4,1,1),oma=c(3,1,1,1))


plot(TotRec[1:(yearsDat-1)],type="l",las=1,ylab="Recruitment")
lines(trueRec[1:(yearsDat-1)],lty=2,col=2)
plot(TotFdir,type="l",las=1,ylab="Directed F")
lines(trueFmort,lty=2,col=2)
barplot(exp(stNatLen))
#==plot estimated selectivities

cnt<-grep("srv_sel50",data2)
sel50surv<-as.numeric(data2[cnt+1])
cnt<-grep("srv_sel95",data2)
sel95surv<-as.numeric(data2[cnt+1])

cnt<-grep("fish_sel50",data2)
sel50fish<-as.numeric(data2[cnt+1])
cnt<-grep("fish_sel95",data2)
sel95fish<-as.numeric(data2[cnt+1])

SurvSel<-1/( 1 + exp( -1*log(19)*(AgeBin-sel50surv)/(sel95surv-sel50surv)))
FishSel<-1/( 1 + exp( -1*log(19)*(AgeBin-sel50fish)/(sel95fish-sel50fish)))

trueSurvSel<-1/( 1 + exp( -1*log(19)*(LengthBin-trueSurvSelPar[1])/(trueSurvSelPar[2]-trueSurvSelPar[1])))
plot(trueSurvSel~LengthBin,type="l",xlab="Age",ylab="Selectivity")
trueFishSel<-1/( 1 + exp( -1*log(19)*(LengthBin-trueFishSelPar[1])/(trueFishSelPar[2]-trueFishSelPar[1])))

lines(FishSel~AgeBin,lty=2,col=2)
lines(SurvSel~AgeBin,lty=2,col=2)
lines(trueFishSel~LengthBin,lty=2,col=1)
legend("bottomright",bty='n',col=c(1,1,2,2),lty=c(1,2,1,2),legend=c("True Survey","Est Survey","True Fishery","Est Fishery"))

dev.new(width=7,height=6)
par(mfrow=c(2,2),mar=c(.1,.1,.1,.1),oma=c(5,5,1,5))
#plot(obsSurvBio,yaxt='n',ylim=c(0,max(obsSurvBio,predSurvBio,trueSpbio,na.rm=T)))
#axis(side=2,las=1)
#lines(predSurvBio/2)
#lines(trueSpbio,col=2,lty=2)
#lines(predSpBio,col=1,lty=1)
#legend("topright",bty='n',"Survey")

#plot(obsSurvBio,yaxt='n',ylim=c(0,max(obsSurvBio,predSurvBio,trueSpbio,na.rm=T)))
plot(obsSurvBio,yaxt='n',xaxt='n')
axis(side=2,las=1)
lines(predSurvBio)
legend("topright",bty='n',"Survey")

#plot(cpueIndex,yaxt='n',ylim=c(0,max(predCpueBio ,cpueIndex)))
plot(cpueIndex,yaxt='n',xaxt='n')
axis(side=4,las=1)
lines(predCpueBio  )
#lines(trueCatch,col=2,lty=2)
legend("topright",bty='n',"CPUE")

plot(obsCatchBio ,yaxt='n')
lines(predCatchBio,yaxt='n')
axis(side=2,las=1)
legend("topright",bty='n',"Catch")

plot(predSpBio,yaxt='n',type="l",ylim=c(0,max(predSpBio,na.rmm=T)))
lines(trueSpbio,col=2,lty=2)
abline(h=B35,lty=3,col=3)
axis(side=4,las=1)
legend("topright",bty='n',c("Spawning biomass","True","Estimated"),lty=c(NA,1,2))

dev.new()
putIn<-ceiling(sqrt(yearsDat))
par(mfrow=c(putIn,putIn),mar=c(0.1,.1,.1,.1))
for(i in 2:yearsDat)
{
plot(obsSurlen[i,],axes=F)
lines(predSurlen[i,])
}
plot.new()
legend("center","Survey")
dev.new()
putIn<-ceiling(sqrt(yearsDat))
par(mfrow=c(putIn,putIn),mar=c(0.1,.1,.1,.1))
for(i in 2:yearsDat)
{
plot(obsCatchlen[i,],axes=F)
lines(predCatchlen[i,])
}
plot.new()
legend("center","Catch")

}
##################################################

TakeOut<-function(SearchTerm,SearchPool)
{
 TakeIndex<-grep(SearchTerm,SearchPool[,3])
 temp<-suppressWarnings(as.numeric(SearchPool[TakeIndex,1]))
 if(is.na(temp))
  temp<-eval(parse(text=SearchPool[TakeIndex,1]))
 return(temp)
}
#input<-"C:/Users/Cody/Desktop/Publications/BlueShark/CTLfile.csv"
ReadCTLfile<-function(input)
{
 OM				<-NULL
 SearchPool			<-as.matrix(read.csv(input))
 OM$Nsim			<-TakeOut("Nsim",SearchPool)
 OM$SimYear			<-TakeOut("SimYear",SearchPool)
 OM$InitYear		<-TakeOut("InitYear",SearchPool)
 OM$AssessmentType	<-TakeOut("AssessmentType",SearchPool)
 OM$FisheryIndepenDat	<-TakeOut("FisheryIndepenDat",SearchPool)
 OM$LifeHistoryPlots	<-TakeOut("LifeHistoryPlots",SearchPool)
 OM$EstimationPlots	<-TakeOut("EstimationPlots",SearchPool)
 OM$AssessmentData	<-TakeOut("AssessmentData",SearchPool)
 OM$CatchCVn		<-TakeOut("CatchCVn",SearchPool)
 OM$CatchCVs		<-TakeOut("CatchCVs",SearchPool)
 OM$CPUEcvN			<-TakeOut("CPUEcvN",SearchPool)
 OM$CPUEcvS			<-TakeOut("CPUEcvS",SearchPool)
 OM$IndexCVn		<-TakeOut("IndexCVn",SearchPool)
 OM$IndexCVs		<-TakeOut("IndexCVs",SearchPool)
 OM$LenSampleN		<-TakeOut("LenSampleN",SearchPool)
 OM$LenSampleS		<-TakeOut("LenSampleS",SearchPool)
 OM$GrowthSDn		<-TakeOut("GrowthSDn",SearchPool)
 OM$GrowthSDs		<-TakeOut("GrowthSDs",SearchPool)

 OM$DiffFlag		<-TakeOut("DiffFlag",SearchPool)

 OM$MaxAge			<-TakeOut("MaxAge",SearchPool)
 OM$NatMn			<-TakeOut("NatMn",SearchPool)
 OM$VonKn			<-TakeOut("VonKn",SearchPool)
 OM$LinfN			<-TakeOut("LinfN",SearchPool)
 OM$t0n			<-TakeOut("t0n",SearchPool)
 OM$mat50n			<-TakeOut("mat50n",SearchPool)
 OM$mat95n			<-TakeOut("mat95n",SearchPool)
 OM$alphaN			<-TakeOut("alphaN",SearchPool)
 OM$betaN			<-TakeOut("betaN",SearchPool)
 OM$MaxMovingN		<-TakeOut("MaxMovingN",SearchPool)
 OM$Move50n			<-TakeOut("Move50n",SearchPool)
 OM$Move95n			<-TakeOut("Move95n",SearchPool)
 OM$steepnessN		<-TakeOut("steepnessN",SearchPool)
 OM$RzeroN			<-TakeOut("RzeroN",SearchPool)
 OM$sigmaRn			<-TakeOut("sigmaRn",SearchPool)

 OM$sel50n			<-TakeOut("sel50n",SearchPool)
 OM$sel95n			<-TakeOut("sel95n",SearchPool)
 OM$surv50n			<-TakeOut("surv50n",SearchPool)
 OM$surv95n			<-TakeOut("surv95n",SearchPool)

 OM$NatMs			<-TakeOut("NatMs",SearchPool)
 OM$VonKs			<-TakeOut("VonKs",SearchPool)
 OM$LinfS			<-TakeOut("LinfS",SearchPool)
 OM$t0s			<-TakeOut("t0s",SearchPool)
 OM$mat50s			<-TakeOut("mat50s",SearchPool)
 OM$mat95s			<-TakeOut("mat95s",SearchPool)
 OM$alphaS			<-TakeOut("alphaS",SearchPool)
 OM$betaS			<-TakeOut("betaS",SearchPool)
 OM$MaxMovingS		<-TakeOut("MaxMovingS",SearchPool)
 OM$Move50s			<-TakeOut("Move50s",SearchPool)
 OM$Move95s			<-TakeOut("Move95s",SearchPool)
 OM$steepnessS		<-TakeOut("steepnessS",SearchPool)
 OM$RzeroS			<-TakeOut("RzeroS",SearchPool)
 OM$sigmaRs			<-TakeOut("sigmaRs",SearchPool)

 OM$sel50s			<-TakeOut("sel50s",SearchPool)
 OM$sel95s			<-TakeOut("sel95s",SearchPool)
 OM$surv50s			<-TakeOut("surv50s",SearchPool)
 OM$surv95s			<-TakeOut("surv95s",SearchPool)

 OM$HistoricalF		<-TakeOut("HistoricalF",SearchPool)
 OM$CatchShareN		<-TakeOut("CatchShareN",SearchPool)
 OM$depletion		<-TakeOut("depletion",SearchPool)
 OM$HarvestControl	<-TakeOut("HarvestControl",SearchPool)
 OM$ConstantCatch		<-TakeOut("ConstantCatch",SearchPool)
 OM$ConstantF		<-TakeOut("ConstantF",SearchPool)
 OM$HCalpha			<-TakeOut("HCalpha",SearchPool)
 OM$HCbeta			<-TakeOut("HCbeta",SearchPool)

 OM$SmalNum			<-TakeOut("SmallNum",SearchPool)
 OM$InitSmooth		<-TakeOut("InitSmooth",SearchPool)
 OM$FmortPen		<-TakeOut("FmortPen",SearchPool)
 OM$RecruitPen		<-TakeOut("RecruitPen",SearchPool)

 list(OM=OM)
}
#====================================================
CleanInput<-function(input,Simyear)
{
if(length(input)==1)
 output		<-rep(input,SimYear)
if(length(input)>1)
 output		<-input
 return(output)
}
#=============================================
PlotLifeHistory<-function()
  {
 dev.new()
 par(mfrow=c(4,4),mar=c(3,3,0,0),oma=c(1,3,1,1))


 plot(LenAtAgeN[1,],type="l",las=1)
 for(x in 2:SimYear)
  lines(LenAtAgeN[x,])
 mtext(side=3,"Length at Age (N)",cex=.7)

 plot(LenAtAgeS[1,],type="l",las=1)
 for(x in 2:SimYear)
  lines(LenAtAgeS[x,])
 mtext(side=3,"Length at Age (S)",cex=.7)

 plot(matureN[1,],las=1,type="l")
 for(x in 2:SimYear)
  lines(matureN[x,])
 mtext(side=3,"Maturity at age (N)",cex=.7)

 plot(matureS[1,],las=1,type="l")
 for(x in 2:SimYear)
  lines(matureS[x,])
 mtext(side=3,"Maturity at age (S)",cex=.7)

 plot(vulnN[1,],las=1,type="l")
 for(x in 2:SimYear)
  lines(vulnN[x,])
 mtext(side=3,"Fishery selectivity at age(N)",cex=.7)

 plot(vulnS[1,],las=1,type="l")
 for(x in 2:SimYear)
  lines(vulnS[x,])
 mtext(side=3,"Fishery selectivity at age (S)",cex=.7)

 plot(survSelN[1,],las=1,type="l")
 for(x in 2:SimYear)
  lines(survSelN[x,])
 mtext(side=3,"Index selectivity at age(N)",cex=.7)

 plot(survSelS[1,],las=1,type="l")
 for(x in 2:SimYear)
  lines(survSelS[x,])
 mtext(side=3,"Index selectivity at age (S)",cex=.7)

 plot(WeightAtAgeN[1,],las=1,type="l")
 for(x in 2:SimYear)
  lines(WeightAtAgeN[x,])
 mtext(side=3,"Weight at age (N)",cex=.7)

 plot(WeightAtAgeS[1,],las=1,type="l")
 for(x in 2:SimYear)
  lines(WeightAtAgeS[x,])
 mtext(side=3,"Weight at age (S)",cex=.7)

 plot(MovementN[1,],las=1,type="l")
 for(x in 2:SimYear)
  lines(MovementN[x,])
  mtext(side=3,"Movement at age (N)",cex=.7)

 plot(MovementS[1,],las=1,type="l")
 for(x in 2:SimYear)
  lines(MovementS[x,])
  mtext(side=3,"Movement at age (S)",cex=.7)

ssb<-seq(1,VirBioN,VirBioN/100)
record<-ssb
  for(x in 1:length(record))
   {
     record[x]<-Recruitment(EggsIN=ssb[x],steepnessIN=steepnessN[1],RzeroIN=RzeroN[1],RecErrIN=RecErrN[1],recType="BH",NatMin=NatMn[1],
							vulnIN=vulnN[1,],matureIN=matureN[1,],weightIN=WeightAtAgeN[1,],LenAtAgeIN=LenAtAgeN[1,])
   }
plot(record~ssb)
mtext(side=3,"SSB vs Rec (N)",cex=.7)
ssb<-seq(1,VirBioS,VirBioS/100)
record<-ssb
  for(x in 1:length(record))
   {
     record[x]<-Recruitment(EggsIN=ssb[x],steepnessIN=steepnessS[1],RzeroIN=RzeroS[1],RecErrIN=RecErrS[1],recType="BH",NatMin=NatMs[1],
							vulnIN=vulnS[1,],matureIN=matureS[1,],weightIN=WeightAtAgeS[1,],LenAtAgeIN=LenAtAgeS[1,])
   }
plot(record~ssb)
mtext(side=3,"SSB vs Rec (S)",cex=.7)


plot(HistoricalF)
 }

#=======================================================

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
    
 
 

