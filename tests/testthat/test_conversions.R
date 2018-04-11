# library(rcosmo)
# library(Rcpp)
# library(testthat)


## NOTE: NEW REFERENCE FILES ARE CREATED AUTOMATICALLY IF THE OLD ONES ARE DELETED
## LAST REFERENCE FILE CREATION DATE: 12/03/2018


# pix2coords --------------------------------------------------------------


testthat::context("Convert HEALPix to spherical coordinates")

testthat::test_that("Output matrix is as expected for 'ring' ordering", {
  ns <- 2^c(1,2,3,4,5)
  filenames <- paste0("references/pix2coords_ring_spherical_", formatC(ns, width=2, flag="0"), ".rds")
  for (i in 1:length(ns))
  {
    eval(bquote(testthat::expect_equal_to_reference(pix2coords_internal(nside = ns[.(i)],
                                                              nested = FALSE,
                                                              cartesian = FALSE),
                                                    filenames[.(i)])))
  }
})

testthat::test_that("Output matrix is as expected for 'nested' ordering", {
  ns <- 2^c(1,2,3,4,5)
  filenames <- paste0("references/pix2coords_nested_spherical_", formatC(ns, width=2, flag="0"), ".rds")
  for (i in 1:length(ns))
  {
    eval(bquote(testthat::expect_equal_to_reference(pix2coords_internal(nside = ns[.(i)],
                                                               nested = TRUE,
                                                               cartesian = FALSE),
                                                   filenames[.(i)])))
  }
})

testthat::test_that("random sample (spix) results agree with full results (spherical) (nested/ring)", {
  for (ns in 2^c(1,2,3,4,5)){
    for (i in sample(12*ns^2, 10)){
      eval(bquote(testthat::expect_equal(pix2coords_internal(nside = .(ns),
                                                    nested = TRUE,
                                                    spix = .(i),
                                                    cartesian = FALSE)[1,],
                                         pix2coords_internal(nside = .(ns),
                                                    nested = TRUE,
                                                    cartesian = FALSE)[.(i),]
                                         )))
      eval(bquote(testthat::expect_equal(pix2coords_internal(nside = .(ns),
                                                    nested = FALSE,
                                                    spix = .(i),
                                                    cartesian = FALSE)[1,],
                                         pix2coords_internal(nside = .(ns),
                                                    nested = FALSE,
                                                    cartesian = FALSE)[.(i),]
                                         )))
    }
  }
})

testthat::context("Convert HEALPix to cartesian coordinates")

testthat::test_that("Output matrix is as expected for 'ring' ordering", {
  ns <- 2^c(1,2,3,4,5)
  filenames <- paste0("references/pix2coords_ring_cartesian_", formatC(ns, width=2, flag="0"), ".rds")
  for (i in 1:length(ns))
  {
    eval(bquote(testthat::expect_equal_to_reference(pix2coords_internal(nside = ns[.(i)],
                                                               nested = FALSE,
                                                               cartesian = TRUE),
                                                     filenames[.(i)])))
  }
})

testthat::test_that("Output matrix is as expected for 'nested' ordering", {
  ns <- 2^c(1,2,3,4,5)
  filenames <- paste0("references/pix2coords_nested_cartesian_", formatC(ns, width=2, flag="0"), ".rds")
  for (i in 1:length(ns))
  {
    eval(bquote(testthat::expect_equal_to_reference(pix2coords_internal(nside = ns[.(i)],
                                                               nested = TRUE,
                                                               cartesian = TRUE),
                                                    filenames[.(i)])))
  }
})

testthat::test_that("random sample (spix) results agree with full results (cartsian) (nested/ring)", {
  for (ns in 2^c(1,2,3,4,5)){
    for (i in sample(12*ns^2, 10)){
      eval(bquote(testthat::expect_equal(pix2coords_internal(nside = .(ns),
                                                    nested = TRUE,
                                                    spix = .(i),
                                                    cartesian = TRUE)[1,],
                                         pix2coords_internal(nside = .(ns),
                                                    nested = TRUE,
                                                    cartesian = TRUE)[.(i),]
      )))
      eval(bquote(testthat::expect_equal(pix2coords_internal(nside = .(ns),
                                                    nested = FALSE,
                                                    spix = .(i),
                                                    cartesian = TRUE)[1,],
                                         pix2coords_internal(nside = .(ns),
                                                    nested = FALSE,
                                                    cartesian = TRUE)[.(i),]
      )))
    }
  }
})


testthat::context("Compare pix2coord vs python's HEALPy result")

## CODE USED TO GENERATE CARTESIAN TEST DATA FROM .mat HEALPy DATA ON 12/03/2018:
#
# for (Nside in 2^c(1,2,3,4,5))
# {
#   # NEST CARTESIAN
#   mn <- paste0("exploration/tests/verify_py/nest_ns",toString(Nside),".mat")
#   Pxn <- readMat(mn)
#   Pxn <- t(Pxn$hpt.nest)
#   saveRDS(Pxn, file = paste0("tests/testthat/testdata/verify_py/cartesian_nest_ns",
#                              formatC(Nside, width=2, flag="0"), ".rds"))
#   # RING CARTESIAN
#   mr <- paste0("exploration/tests/verify_py/ring_ns",toString(Nside),".mat")
#   Pxr <- readMat(mr)
#   Pxr <- t(Pxr$hpt.ring)
#   saveRDS(Pxr, file = paste0("tests/testthat/testdata/verify_py/cartesian_ring_ns",
#                              formatC(Nside, width=2, flag="0"), ".rds") )
# }

testthat::test_that("cartesian 'ring' data agrees with HEALpy result", {
  ns <- 2^c(1,2,3,4,5)
  filenames <- paste0("testdata/verify_py/cartesian_ring_ns",
                      formatC(ns, width=2, flag="0"), ".rds")
  for (i in 1:length(ns))
  {
    eval(bquote(testthat::expect_equal( readRDS(filenames[.(i)]),
                                        pix2coords_internal(nside = ns[.(i)],
                                                   nested = FALSE,
                                                   cartesian = TRUE)[,1:3] )))
  }
})


testthat::test_that("cartesian 'nested' data agrees with HEALpy result", {
  ns <- 2^c(1,2,3,4,5)
  filenames <- paste0("testdata/verify_py/cartesian_nest_ns",
                      formatC(ns, width=2, flag="0"), ".rds")
  for (i in 1:length(ns))
  {
    eval(bquote(testthat::expect_equal( readRDS(filenames[.(i)]),
                                        pix2coords_internal(nside = ns[.(i)],
                                                   nested = TRUE,
                                                   cartesian = TRUE)[,1:3] )))
  }
})


testthat::test_that("spherical 'ring' data agrees with python results after conversion to cartesian", {
  ns <- 2^c(1,2,3,4,5)
  filenames <- paste0("testdata/verify_py/cartesian_ring_ns",
                      formatC(ns, width=2, flag="0"), ".rds")
  for (i in 1:length(ns))
  {
    p2c <- pix2coords_internal(nside = ns[i], nested = FALSE, cartesian = FALSE)[,1:2]
    p2c <- as.data.frame(p2c)
    names(p2c) <- c("theta", "phi")
    p2c <- as.matrix(sph2car(p2c))
    eval(bquote(testthat::expect_equal( readRDS(filenames[.(i)]), p2c,
                                        check.attributes = FALSE )))
  }
})


testthat::test_that("spherical 'nested' data agrees with python results after conversion to cartesian", {
  ns <- 2^c(1,2,3,4,5)
  filenames <- paste0("testdata/verify_py/cartesian_nest_ns",
                      formatC(ns, width=2, flag="0"), ".rds")
  for (i in 1:length(ns))
  {
    p2c <- pix2coords_internal(nside = ns[i], nested = TRUE, cartesian = FALSE)[,1:2]
    p2c <- as.data.frame(p2c)
    names(p2c) <- c("theta", "phi")
    p2c <- as.matrix(sph2car(p2c))
    eval(bquote(testthat::expect_equal( readRDS(filenames[.(i)]), p2c,
                                        check.attributes = FALSE )))
  }
})


