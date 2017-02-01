rm(list=ls())
CurDir<-"C:/Users/Cody/Desktop/GeMS/"
setwd(CurDir)
source("General_MSE_helpers.R")
source("General_MSE_main.R")

CreateFolderNameList<-c("BlueShark_Base","BlueShark_Move")
CreateFolderNameList<-c("BlueShark_Move")
				
#==Loop that reads in CTL files stored in CurDir and executes code
for(x in 1:length(CreateFolderNameList))
{
 Inout<-ReadCTLfile(paste(CurDir,CreateFolderNameList[x],"_CTL.csv",sep=""))
 GeMS(out=Inout,CreateFolderName=CreateFolderNameList[x])
 ProductionModelOutput(CreateFolderNameList[x])
}


