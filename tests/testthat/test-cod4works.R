context("Base model in Cod_4_MPA runs")

ExampleFilePath<-system.file("extdata",package="GeMS")
truedir<-file.path(ExampleFilePath, "Cod_4_MPA")
MSEdir<-system.file("tests",package="GeMS")

CTLfilename<-"Cod_MPA_CTL"
CTLfilepath<-file.path(truedir,paste0(CTLfilename,".csv"))

test_that("Cod_MPA exists.", {
  expect_true(file.copy(from=CTLfilepath,to=MSEdir,overwrite=T))
})

CTLfile<-ReadCTLfile(file.path(MSEdir,CTLfilename))

test_that("Cod_MPA_CTL reads in.", {
  expect_equal(sum(is.na(CTLfile$OM)), 0)
})

test_that("Cod_MPA_CTL runs.", {
  newCTLfile<-CTLfile
  newCTLfile$OM$Nsim<-1
  newCTLfile$OM$SimYear<-newCTLfile$OM$InitYear+1
  GeMS(newCTLfile,CTLfilename,MSEdir,silent=T)
  Dirname<- file.path(MSEdir,CTLfilename)
  expect_equal(sum(file.exists(file.path(Dirname,"ProdOutputs.csv"))),1)
  unlink(file.path(MSEdir,paste0(CTLfilename,".csv")))
  unlink(Dirname,recursive=T)
  unlink(file.path(MSEdir,"plots"),recursive=T)
})
