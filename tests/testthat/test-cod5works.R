context("Base model in Cod_5_AgeStructure runs")

ExampleFilePath<-system.file("extdata",package="GeMS")
truedir<-file.path(ExampleFilePath, "Cod_5_AgeStructure")
MSEdir<-system.file("tests",package="GeMS")

CTLfilename<-"Cod_AgeStructure_CTL"
CTLfilepath<-file.path(truedir,paste0(CTLfilename,".csv"))

test_that("Cod_AgeStructure_CTL exists.", {
  expect_true(file.copy(from=CTLfilepath,to=MSEdir,overwrite=T))
})

CTLfile<-ReadCTLfile(file.path(MSEdir,CTLfilename))

test_that("Cod_AgeStructure_CTL reads in.", {
  expect_equal(sum(is.na(CTLfile$OM)), 0)
})

test_that("Cod_AgeStructure_CTL runs.", {
  newCTLfile<-CTLfile
  newCTLfile$OM$Nsim<-1
  newCTLfile$OM$SimYear<-newCTLfile$OM$InitYear+1
  GeMS(newCTLfile,CTLfilename,MSEdir,silent=T)
  Dirname<- file.path(MSEdir,CTLfilename)
  expect_equal(sum(file.exists(file.path(Dirname,newCTLfile$OM$Nsim,newCTLfile$OM$SimYear,"SimAss.rep"))),1)
  unlink(file.path(MSEdir,paste0(CTLfilename,".csv")))
  unlink(file.path(Dirname),recursive=T)
  unlink(file.path(MSEdir,"plots"),recursive=T)
})
