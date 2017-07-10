# library(rcosmo)
# library(Rcpp)
# library(testthat)

testthat::context("Convert HEALPix to spherical coordinates and (j,i) indices")

testthat::test_that("Output matrix is as expected for RING ordering with small Nside", {
  testthat::expect_equal_to_reference(pix2angC(2,FALSE), "references/pix2angRING_02.rds")
  testthat::expect_equal_to_reference(pix2angC(4,FALSE), "references/pix2angRING_04.rds")
  testthat::expect_equal_to_reference(pix2angC(8,FALSE), "references/pix2angRING_08.rds")
  testthat::expect_equal_to_reference(pix2angC(16,FALSE), "references/pix2angRING_16.rds")
  testthat::expect_equal_to_reference(pix2angC(32,FALSE), "references/pix2angRING_32.rds")
})

testthat::test_that("Output matrix is as expected for NEST ordering with small Nside", {
  testthat::expect_equal_to_reference(pix2angC(2,TRUE), "references/pix2angNEST_02.rds")
  testthat::expect_equal_to_reference(pix2angC(4,TRUE), "references/pix2angNEST_04.rds")
  testthat::expect_equal_to_reference(pix2angC(8,TRUE), "references/pix2angNEST_08.rds")
  testthat::expect_equal_to_reference(pix2angC(16,TRUE), "references/pix2angNEST_16.rds")
  testthat::expect_equal_to_reference(pix2angC(32,TRUE), "references/pix2angNEST_32.rds")
})

testthat::test_that("Sample pixel results agree with full pixel results with small Nside", {
  for (Ns in c(2,4,8,16,32)){
    for (i in c(1,20,12*Ns^2)){
      eval(bquote(testthat::expect_equal(pix2angC(.(Ns),TRUE,.(i))[1,],pix2angC(.(Ns),TRUE)[.(i),])))
      eval(bquote(testthat::expect_equal(pix2angC(.(Ns),FALSE,.(i))[1,],pix2angC(.(Ns),FALSE)[.(i),])))
    }
  }
})
