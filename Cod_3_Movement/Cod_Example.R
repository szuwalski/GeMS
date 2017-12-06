rm(list=ls())
GeMS.dir <- "C:/GeMS"
Cur.dir<-"C:/GeMS/Cod_3_Movement"

source(file.path(GeMS.dir,"run_GeMS.R"))

OMNames<-c("Cod_Movement_CTL")

run_GeMS(GeMSDir=GeMS.dir, CurDir=Cur.dir,
         CreateFolderNameList=OMNames)