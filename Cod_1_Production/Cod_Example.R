rm(list=ls())

GeMS.dir <- "C:/GeMS"
Cur.dir<-"C:/GeMS/Cod_1_Production/"

source(file.path(GeMS.dir,"run_GeMS.R"))

OMNames<-c("Cod_Base_CTL","Cod_LowProd_CTL")

run_GeMS(GeMSDir=GeMS.dir, CurDir=Cur.dir,
         CreateFolderNameList=OMNames)

##-----------------
## Parallel example
##-----------------
#run_GeMS(GeMSDir=GeMS.dir, CurDir=Cur.dir,
#         CreateFolderNameList=OMNames,
#         runparallel = T, cores = 2)
