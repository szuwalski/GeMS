library(GeMS)
Cur.dir<-system.file("extdata/Cod_4_MPA",package="GeMS")
OMNames<-c("Cod_MPA_CTL")

run_GeMS(MSEdir=Cur.dir,
         CTLNameList=OMNames)
