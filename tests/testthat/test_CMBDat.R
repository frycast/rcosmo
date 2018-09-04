testthat::context("Read CMB data from a FITS file")

testthat::test_that("Output is expected with small example CMB map with Nside 1024 and 5 columns", {
  testthat::expect_equal_to_reference(CMBDat("testdata/CMB_testmap_1024_5cols_10rows.fits", spix = 1:10),
                                      "references/readFITScmb_1204_5cols_10rows.rds")
})
