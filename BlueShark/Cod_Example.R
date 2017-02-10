rm(list=ls())
source("C:/Users/liweiwen/Desktop/General-MSE-1/General-MSE-1-master/General_MSE_helpers.R")
source("C:/Users/liweiwen/Desktop/General-MSE-1/General-MSE-1-master/General_MSE_main.R")

CurDir<-"C:/Users/liweiwen/Desktop/General-MSE-1/General-MSE-1-master/Cod_1_Production/"
setwd(CurDir)
CreateFolderNameList<-c("Cod_Base_CTL","Cod_LowProd_CTL")

#==Loop that reads in CTL files stored in CurDir and executes code
for(x in 1:length(CreateFolderNameList))
{
 Inout<-ReadCTLfile(paste(CurDir,CreateFolderNameList[x],".csv",sep=""))
 GeMS(out=Inout,CreateFolderName=CreateFolderNameList[x])
 ProductionModelOutput(CreateFolderNameList[x])
}



Inout<-ReadCTLfile(paste(CurDir,CreateFolderNameList[1],".csv",sep=""))

#==read in production results==
Data<-readLines(paste(CurDir,CreateFolderNameList[1],"/ProdOutputs.csv",sep=""))

#==estimated and true CPUE
tmp<-grep("est cpue",Data)
Quant<-Data[(tmp+1):(tmp+Inout$OM$Nsim)]
estCPUE<-matrix(ncol=Inout$OM$SimYear,nrow=Inout$OM$Nsim)

for(x in 1:length(Quant))
 estCPUE[x,]<-as.numeric(unlist(strsplit(Quant[x],split=" ")))

tmp<-grep("true CPUE",Data)
Quant<-Data[(tmp+1):(tmp+Inout$OM$Nsim)]
trueCPUE<-matrix(ncol=Inout$OM$SimYear,nrow=Inout$OM$Nsim)
for(x in 1:length(Quant))
 trueCPUE[x,]<-as.numeric(unlist(strsplit(Quant[x],split=" ")))

plot(trueCPUE[1,],type="l",ylim=c(0,max(trueCPUE,na.rm=T)),las=1,yaxt='n',ylab='')
for(x in 1:nrow(estCPUE))
 lines(estCPUE[x,])
lines(trueCPUE[1,],col=2,lwd=2)

#==True BMSY==
tmp<-grep("true BMSY",Data)
trueBMSY<-as.numeric(Data[(tmp+1)])

tmp<-grep("est BMSY",Data)
estBMSY<-matrix(ncol=Inout$OM$SimYear,nrow=Inout$OM$Nsim)
for(x in 1:length(Quant))
 estBMSY[x,]<-as.numeric(unlist(strsplit(Quant[x],split=" ")))

temp<-sweep(estBMSY,MAR=2,trueBMSY,FUN="-")
RelativeErrorBMSY<-sweep(temp,MAR=2,trueBMSY,FUN="/")
boxplot(RelativeErrorBMSY[,Inout$OM$InitYear:Inout$OM$SimYear],las=1)
mtext(side=2,"Relative error",line=2.5)

#==estimated and true biomass###
tmp<-grep("est spawning biomass",Data)
Quant<-Data[(tmp+1):(tmp+Inout$OM$Nsim)]
estBIOMA<-matrix(ncol=Inout$OM$SimYear,nrow=Inout$OM$Nsim)

for(x in 1:length(Quant))
 estBIOMA[x,]<-as.numeric(unlist(strsplit(Quant[x],split=" ")))

tmp<-grep("true spawning biomass",Data)
Quant<-Data[(tmp+1):(tmp+Inout$OM$Nsim)]
trueBIOMA<-matrix(ncol=Inout$OM$SimYear,nrow=Inout$OM$Nsim)
for(x in 1:length(Quant))
 trueBIOMA[x,]<-as.numeric(unlist(strsplit(Quant[x],split=" ")))

plot(trueBIOMA[1,],type="l",ylim=c(0,max(estBIOMA,na.rm=T)),las=1,yaxt='n',ylab='')
for(x in 1:nrow(estBIOMA))
 points(estBIOMA[x,],pch=16)
lines(trueBIOMA[1,],col=2,lwd=3)

###relative error of BMSY###
tmp<-grep("true BMSY",Data)
trueBMSY<-as.numeric(Data[(tmp+1)])

tmp<-grep("est BMSY",Data)
estBMSY<-matrix(ncol=Inout$OM$SimYear,nrow=Inout$OM$Nsim)
Quant<-Data[(tmp+1):(tmp+Inout$OM$Nsim)]
for(x in 1:length(Quant))
 estBMSY[x,]<-as.numeric(unlist(strsplit(Quant[x],split=" ")))
relaerror<-(estBMSY-trueBMSY)/trueBMSY
boxplot(relaerror)

###relative error of FMSY###

tmp<-grep("true FMSY",Data)
trueFMSY<-as.numeric(Data[(tmp+1)])

tmp<-grep("est FMSY",Data)
estFMSY<-matrix(ncol=Inout$OM$SimYear,nrow=Inout$OM$Nsim)
Quant<-Data[(tmp+1):(tmp+Inout$OM$Nsim)]
for(x in 1:length(Quant))
 estFMSY[x,]<-as.numeric(unlist(strsplit(Quant[x],split=" ")))
relaerror<-(estFMSY-trueFMSY)/trueFMSY
plot(relaerror[1,],type="l",ylim=c(0,max(trueFMSY,na.rm=T)),las=1,yaxt='n',ylab='')
boxplot(relaerror)


###relative error TAC###
tmp<-grep("true total allowable catch",Data)
trueTAC<-as.numeric(Data[(tmp+1)])

tmp<-grep("est total allowable catch",Data)
estTAC<-matrix(ncol=Inout$OM$SimYear,nrow=Inout$OM$Nsim)
Quant<-Data[(tmp+1):(tmp+Inout$OM$Nsim)]
for(x in 1:length(Quant))
 estTAC[x,]<-as.numeric(unlist(strsplit(Quant[x],split=" ")))
relaerror<-(estTAC-trueTAC)/trueTAC   
plot(relaerror[1,],type="l",ylim=c(0,max(trueTAC,na.rm=T)),las=1,yaxt='n',ylab='')

###Relative error in estimated terminal year biomass over time###

tmp<-grep("true spawning biomass",Data)
trueSPW<-as.numeric(Data[(tmp+1)])

tmp<-grep("est spawning biomass",Data)
estSPW<-matrix(ncol=Inout$OM$SimYear,nrow=Inout$OM$Nsim)
Quant<-Data[(tmp+1):(tmp+Inout$OM$Nsim)]
for(x in 1:length(Quant))
 estSPW[x,]<-as.numeric(unlist(strsplit(Quant[x],split=" ")))
relaerror<-(estSPW-trueSPW)/trueSPW   
plot(relaerror[1,],type="l",ylim=c(0,max(trueSPW,na.rm=T)),las=1,yaxt='n',ylab='')


###¡®true¡¯ stock status###

tmp<-grep("true spawning biomass",Data)
Quant<-Data[(tmp+1):(tmp+Inout$OM$Nsim)]
trueBIOMA<-matrix(ncol=Inout$OM$SimYear,nrow=Inout$OM$Nsim)
for(x in 1:length(Quant))
 trueBIOMA[x,]<-as.numeric(unlist(strsplit(Quant[x],split=" ")))

tmp<-grep("true BMSY",Data)
trueBMSY<-as.numeric(Data[(tmp+1)])

status<-trueBIOMA/trueBMSY
boxplot(status,type="1",las=1,ylab='',ylim=c(0,max(status,na.rm=T)))
abline(h=1,lty=2,col=2)

###Probability of an estimated overfished stock vs. the probability of a true overfished stock###
tmp<-grep("true BMSY",Data)
trueBMSY<-as.numeric(Data[(tmp+1)])

tmp<-grep("est spawning biomass",Data)
estSPW<-matrix(ncol=Inout$OM$SimYear,nrow=Inout$OM$Nsim)
Quant<-Data[(tmp+1):(tmp+Inout$OM$Nsim)]
for(x in 1:length(Quant))
 estSPW[x,]<-as.numeric(unlist(strsplit(Quant[x],split=" ")))

p1<-estSPW/trueBMSY
par(new=T)
boxplot(p1,type="1",las=1,ylab='',col=2,ylim=c(0,3))
abline(h=1,lty=2,col=2)

boxplot(status[,50:80],type="1",las=1,ylab='',ylim=c(0,2))
par(new=T)
boxplot(p1[,50:80],type="1",las=1,ylab='',col=2,ylim=c(0,1))
abline(h=1,lty=2,col=2)
