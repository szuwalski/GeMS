context("Base model in Cod_3_Movement runs")

ExampleFilePath<-system.file("extdata",package="GeMS")
MSEdir<-file.path(ExampleFilePath, "Cod_3_Movement")
CTLfilename<-"Cod_Movement_CTL"
CTLfile<-ReadCTLfile(file.path(MSEdir,CTLfilename))

test_that("Cod_Movement_CTL reads in.", {
  expect_equal(sum(is.na(CTLfile$OM)), 0)
})

test_that("Cod_Movement_CTL runs.", {
  newCTLfile<-CTLfile
  newCTLfile$OM$Nsim<-1
  newCTLfile$OM$SimYear<-newCTLfile$OM$InitYear+1
  GeMS(newCTLfile,CTLfilename,MSEdir)
  PNGname<- file.path(MSEdir,"plots",paste0("ProductionRefPoints_",CTLfilename,".png"))
  expect_equal(sum(file.exists(PNGname)),1)
  unlink(file.path(MSEdir,CTLfilename),recursive=T)
  unlink(file.path(MSEdir,"plots"),recursive=T)
})
