rm(list=ls())
Cur.dir<-system.file("Cod_3_Movement",package="GeMS")

OMNames<-c("Cod_Movement_CTL")

run_GeMS(MSEdir=Cur.dir,
         CTLNameList=OMNames)