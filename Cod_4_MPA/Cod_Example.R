rm(list=ls())
GeMS.dir <- "C:/GeMS"
Cur.dir<-"C:/GeMS/Cod_4_MPA/"

source(file.path(GeMS.dir,"run_GeMS.R"))
OMNames<-c("Cod_MPA_CTL")

run_GeMS(GeMSDir=GeMS.dir, CurDir=Cur.dir,
         CreateFolderNameList=OMNames)
