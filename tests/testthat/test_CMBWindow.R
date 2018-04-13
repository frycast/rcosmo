testthat::context("area of CMBDataFrame and CMBWindow")

# Create windows and corresponding cmbdfs of each winType and coords
cmbdf <- CMBDataFrame(nside = 128, ordering = "nested",
                      coords = "cartesian")

disc.xyz <- CMBWindow(x = 1, y = 0, z = 0, r = 0.5)
m.disc.xyz <- CMBWindow(x = 1, y = 0, z = 0, r = 0.5, set.minus = TRUE)
polygon <- CMBWindow(phi = c(0, pi/4, pi/4, pi/5),
                     theta = c(pi/2, pi/2, pi/4, pi/2 - pi/20))
m.polygon <- CMBWindow(phi = c(0, pi/4, pi/4, pi/5),
                       theta = c(pi/2, pi/2, pi/4, pi/2 - pi/20),
                       set.minus = TRUE)

disc <- coords(disc.xyz, new.coords = "spherical")
m.disc <- coords(m.disc.xyz, new.coords = "spherical")
poly.xyz <- coords(polygon, new.coords = "cartesian")
m.poly.xyz <- coords(m.polygon, new.coords = "cartesian")

disc.cmb <- window(cmbdf, new.window = disc)
m.disc.cmb <- window(cmbdf, new.window = m.disc)
poly.cmb <- window(cmbdf, new.window = polygon)
m.poly.cmb <- window(cmbdf, new.window = m.polygon)


testthat::test_that("areas sum to 4*pi", {
  testthat::expect_equal(geoArea(disc.cmb) + geoArea(m.disc.cmb), 4*pi)
  testthat::expect_equal(geoArea(poly.cmb) + geoArea(m.poly.cmb), 4*pi)
  testthat::expect_equal(geoArea(disc) + geoArea(m.disc), 4*pi)
  testthat::expect_equal(geoArea(polygon) + geoArea(m.polygon), 4*pi)
})


testthat::test_that("areas differences are accounted for", {
  testthat::expect_equal(  geoArea(poly.cmb) - geoArea(polygon),
                         -(geoArea(m.poly.cmb) - geoArea(m.polygon)))
  testthat::expect_equal(  geoArea(disc.cmb) - geoArea(disc),
                         -(geoArea(m.disc.cmb) - geoArea(m.disc)) )
})



testthat::context("CMBWindows have correct attributes")

testthat::test_that("CMBWindow column names are correct", {

  testthat::expect_equal( names(polygon), c("theta","phi") )
  testthat::expect_equal( names(m.polygon), c("theta","phi") )
  testthat::expect_equal( names(disc), c("theta","phi","r") )
  testthat::expect_equal( names(m.disc), c("theta","phi","r") )

  testthat::expect_equal( names(poly.xyz), c("x","y","z") )
  testthat::expect_equal( names(m.poly.xyz), c("x","y","z") )
  testthat::expect_equal( names(disc.xyz), c("x","y","z","r") )
  testthat::expect_equal( names(m.disc.xyz), c("x","y","z","r") )
})

testthat::test_that("CMBWindow coords names are correct", {

  testthat::expect_equal( coords(polygon), "spherical" )
  testthat::expect_equal( coords(m.polygon), "spherical" )
  testthat::expect_equal( coords(disc), "spherical" )
  testthat::expect_equal( coords(m.disc), "spherical" )

  testthat::expect_equal( coords(poly.xyz), "cartesian" )
  testthat::expect_equal( coords(m.poly.xyz), "cartesian" )
  testthat::expect_equal( coords(disc.xyz), "cartesian" )
  testthat::expect_equal( coords(m.disc.xyz), "cartesian" )
})

testthat::test_that("CMBWindow winTypes are correct", {

  testthat::expect_equal( winType(polygon), "polygon" )
  testthat::expect_equal( winType(m.polygon), "minus.polygon" )
  testthat::expect_equal( winType(disc), "disc" )
  testthat::expect_equal( winType(m.disc), "minus.disc" )

  testthat::expect_equal( winType(poly.xyz), "polygon" )
  testthat::expect_equal( winType(m.poly.xyz), "minus.polygon" )
  testthat::expect_equal( winType(disc.xyz), "disc" )
  testthat::expect_equal( winType(m.disc.xyz), "minus.disc" )

  testthat::expect_equal( winType(list(polygon, m.polygon,
                                       disc, m.disc)),
                          c("polygon", "minus.polygon",
                            "disc", "minus.disc"))

  testthat::expect_equal( winType(
    winType(polygon, new.type = "minus.polygon")
  ), "minus.polygon" )
})

testthat::test_that("CMBWindow assumedConvex is correct", {

  testthat::expect_equal( assumedConvex(polygon), FALSE )
  testthat::expect_equal( assumedConvex(m.polygon), FALSE )
  testthat::expect_equal( assumedConvex(disc), TRUE )
  testthat::expect_equal( assumedConvex(m.disc), TRUE )

  testthat::expect_equal( assumedConvex(poly.xyz), FALSE )
  testthat::expect_equal( assumedConvex(m.poly.xyz), FALSE )
  testthat::expect_equal( assumedConvex(disc.xyz), TRUE )
  testthat::expect_equal( assumedConvex(m.disc.xyz), TRUE )

  testthat::expect_equal( assumedConvex(
    assumedConvex(polygon, assume.convex = TRUE)),
    TRUE )

  testthat::expect_equal( assumedConvex(
    CMBWindow(phi = c(0, pi/4, pi/4, pi/5),
              theta = c(pi/2, pi/2, pi/4, pi/2 - pi/20),
              assume.convex = TRUE)),
    TRUE )

})


### MORE TO DO IN FUTURE:
# testthat::context("CMBWindow helper functions")
# e.g. maxDist, triangulate

