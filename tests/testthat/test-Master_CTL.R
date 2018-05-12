context("Master_CTL.csv is keeping up ReadCTLfile()")

test_that("ReadCTL() works with Master_CTL", {
  MasterFilePath<-system.file("extdata",package="GeMS")
  expect_equal(sum(is.na(ReadCTLfile(file.path(MasterFilePath,"Master_CTL"))$OM)), 0)
})