#' Function to plot age-structured models
#' Used within the \code{\link{GeMS}} function, not to be run on its own.
#' 
#' @param REP Read SimAss.rep file
#' @param CTL Read SimAss.CTL file
#' @param DAT Read SimAss.DAT file
#' @param TRU Read TrueQuantities.DAT file
#' @param data2 Read Simass.PAR file
#' @param MSEdir Directory containing CTL files
#' @param CTLNameList Vector of CTL file names
#' 
#' @return GeMS Object # To be implemented
#'
#' @export
#' 
AgeAssPlot<-function(REP,CTL,DAT,TRU,data2,MSEdir,CTLNameList)
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
obsSurlen<-matrix(suppressWarnings(as.numeric(unlist(strsplit(REP[(temp+1):(temp+yearsDat)],split=" ")))),nrow=yearsDat,byrow=T)[,2:(length(LengthBin)+1)]
temp<-grep("pred surv len freq",REP)
predSurlen<-matrix(suppressWarnings(as.numeric(unlist(strsplit(REP[(temp+1):(temp+yearsDat)],split=" ")))),nrow=yearsDat,byrow=T)[,2:(length(LengthBin)+1)]

temp<-grep("obs catch len freq",REP)
obsCatchlen<-matrix(suppressWarnings(as.numeric(unlist(strsplit(REP[(temp+1):(temp+yearsDat)],split=" ")))),nrow=yearsDat,byrow=T)[,2:(length(LengthBin)+1)]
temp<-grep("pred catch len freq",REP)
predCatchlen<-matrix(suppressWarnings(as.numeric(unlist(strsplit(REP[(temp+1):(temp+yearsDat)],split=" ")))),nrow=yearsDat,byrow=T)[,2:(length(LengthBin)+1)]

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



png(file.path(MSEdir,"plots",paste0("PopulationProcessesEst_",paste(CTLNameList,sep="_"),".png")),res=1200,width=6,height=6,units="in")
par(mfrow=c(3,1),mar=c(1,4,1,1),oma=c(3,1,1,1))
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

cnt<-grep("log_avg_SelPars50",data2)
sel50fish<-exp(as.numeric(data2[cnt+1]))
cnt<-grep("log_avg_SelPars95",data2)
sel95fish<-exp(as.numeric(data2[cnt+1]))

cnt<-grep("SelPars_dev50",data2)
sel50_dir_dev<-as.numeric(data2[(cnt+1):(cnt+length(fmortYrs))])
cnt<-grep("SelPars_dev95",data2)
sel95_dir_dev<-as.numeric(data2[(cnt+1):(cnt+length(fmortYrs))])

SurvSel<-1/( 1 + exp( -1*log(19)*(AgeBin-sel50surv)/(sel95surv-sel50surv)))
FishSel<-1/( 1 + exp( -1*log(19)*(AgeBin-sel50fish)/(sel95fish-sel50fish)))

trueSurvSel<-1/( 1 + exp( -1*log(19)*(LengthBin-trueSurvSelPar[1])/(trueSurvSelPar[2]-trueSurvSelPar[1])))
#plot(trueSurvSel~LengthBin,type="l",xlab="Age",ylab="Selectivity")
trueFishSel<-1/( 1 + exp( -1*log(19)*(LengthBin-trueFishSelPar[1])/(trueFishSelPar[2]-trueFishSelPar[1])))
#plot(FishSel~AgeBin,type="l",xlab="Age",ylab="Selectivity")

#lines(FishSel~AgeBin,lty=2,col=2)
#lines(SurvSel~AgeBin,lty=2,col=2)
#lines(trueFishSel~LengthBin,lty=2,col=1)
#legend("bottomright",bty='n',col=c(1,1,2,2),lty=c(1,2,1,2),legend=c("True Survey","Est Survey","True Fishery","Est Fishery"))
dev.off()


png(file.path(MSEdir,"plots", paste0("FitsToDataAgeStruc_",paste(CTLNameList,sep="_",collapse=""),".png")),res=1200,width=6,height=6,units="in")

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
dev.off()



png(file.path(MSEdir,"plots",paste0("FitsToSurvLengths_",paste(CTLNameList,sep="_",collapse=""),".png")),res=1200,width=6,height=6,units="in")
putIn<-ceiling(sqrt(yearsDat))
par(mfrow=c(putIn,putIn),mar=c(0.1,.1,.1,.1))
for(i in 2:yearsDat)
{
plot(obsSurlen[i,],axes=F)
lines(predSurlen[i,])
}
plot.new()
legend("center","Survey")
dev.off()



png(file.path(MSEdir,"plots",paste0("FitsToCatchLengths_",paste(CTLNameList,sep="_",collapse=""),".png")),res=1200,width=6,height=6,units="in")

putIn<-ceiling(sqrt(yearsDat))
par(mfrow=c(putIn,putIn),mar=c(0.1,.1,.1,.1))
for(i in 2:yearsDat)
{
plot(obsCatchlen[i,],axes=F)
lines(predCatchlen[i,])
}
plot.new()
legend("center","Catch")
dev.off()
}