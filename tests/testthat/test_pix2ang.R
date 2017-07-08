library(rcosmo)
library(Rcpp)
library(testthat)
#sourceCpp("src/pix2ang.cpp")

context("Convert HEALPix to spherical coordinates and (j,i) indices")

test_that("Output matrix is as expected for RING ordering with small Nside", {
  expect_equal_to_reference(pix2angC(2,FALSE), "pix2angRING_02.rds")
  expect_equal_to_reference(pix2angC(4,FALSE), "pix2angRING_04.rds")
  expect_equal_to_reference(pix2angC(8,FALSE), "pix2angRING_08.rds")
  expect_equal_to_reference(pix2angC(16,FALSE), "pix2angRING_16.rds")
  expect_equal_to_reference(pix2angC(32,FALSE), "pix2angRING_32.rds")
})

test_that("Output matrix is as expected for NEST ordering with small Nside", {
  expect_equal_to_reference(pix2angC(2,TRUE), "pix2angNEST_02.rds")
  expect_equal_to_reference(pix2angC(4,TRUE), "pix2angNEST_04.rds")
  expect_equal_to_reference(pix2angC(8,TRUE), "pix2angNEST_08.rds")
  expect_equal_to_reference(pix2angC(16,TRUE), "pix2angNEST_16.rds")
  expect_equal_to_reference(pix2angC(32,TRUE), "pix2angNEST_32.rds")
  expect_equal(1,0)
})
