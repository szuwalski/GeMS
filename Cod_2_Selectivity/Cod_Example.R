rm(list=ls())
GeMS.dir <- "C:/GeMS"
Cur.dir<-"C:/GeMS/Cod_2_Selectivity/"

source(file.path(GeMS.dir,"run_GeMS.R"))
OMNames<-c("Cod_Selectivity_CTL")

run_GeMS(GeMSDir=GeMS.dir, CurDir=Cur.dir,
         CreateFolderNameList=OMNames)


