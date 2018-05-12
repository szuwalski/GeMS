library(GeMS)
Cur.dir<-system.file("extdata/Cod_3_Movement",package="GeMS")

OMNames<-c("Cod_Movement_CTL")

run_GeMS(MSEdir=Cur.dir,
         CTLNameList=OMNames)