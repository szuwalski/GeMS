rm(list=ls())
source("C:/GeMS/General_MSE_helpers.R")
source("C:/GeMS/General_MSE_main.R")

CurDir<-"C:/GeMS/BlueShark/"
setwd(CurDir)
CreateFolderNameList<-c("BlueShark_Base")

#==Loop that reads in CTL files stored in CurDir and executes code
for(x in 1:length(CreateFolderNameList))
{
 Inout<-ReadCTLfile(paste(CurDir,CreateFolderNameList[x],".csv",sep=""))
 GeMS(out=Inout,CreateFolderName=CreateFolderNameList[x])
 ProductionModelOutput(CreateFolderNameList[x])
}


