rm(list=ls())
source("C:/GeMS/General_MSE_helpers.R")
source("C:/GeMS/General_MSE_main.R")

CurDir<-"C:/GeMS/Cod_2_Selectivity"
setwd(CurDir)
CreateFolderNameList<-c("Cod_Selectivity_CTL")

#==Loop that reads in CTL files stored in CurDir and executes code
for(x in 1:length(CreateFolderNameList))
{
 Inout<-ReadCTLfile(file.path(CurDir,paste0(CreateFolderNameList[x],".csv")))
 GeMS(out=Inout,CreateFolderName=CreateFolderNameList[x])
 ProductionModelOutput(CreateFolderNameList=CreateFolderNameList[x])
}


