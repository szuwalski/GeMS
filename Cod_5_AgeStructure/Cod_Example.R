rm(list=ls())

#==name where the GeMS folder lives on your computer
dir.GeMS<-"C:/GeMS"

#==source important scripts
source(file.path(TopDir,"run_GeMS.R"))

#==identify the directory where current analysis will be performed and CTL files
dir.MSE<-"C:/GeMS/Cod_5_AgeStructure"

CreateFolderNameList<-c("Cod_AgeStructure_CTL","Cod_Age_Mvary_CTL","Cod_Age_Mvary_estM_CTL")

#==Loop that reads in CTL files stored in dir.MSE and executes code
run_GeMS(GeMSDir=dir.GeMS, CurDir=dir.MSE,
		 CreateFolderNameList=OMNames)

##==compare just two scenarios: estimating M vs. not estimating M when M is varying
##==WARNING: Overwrites plots made from run_GeMS
#CreateFolderNameList<-c("Cod_Age_Mvary_CTL","Cod_Age_Mvary_estM_CTL")
#AgeStructureComp(Inout=Inout,CreateFolderNameList=CreateFolderNameList,MSEdir=dir.MSE)
