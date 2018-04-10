## TESTS NOT YET IMPLEMENTED


# construct a CMBDataFrame, check it has the right class and other attributes
# do this for all types of CMBDataFrame including no coords and NA intensities.
# Then get a window of that CMBDataFrame
# and check that the window has the right attributes etc. Do this for all
# types of windows. Test all the attribute functions here too.
# Do everything with and without unexpected extra columns attached.
# Check that attributes aren't dropped with cmbdf[i,], cmbdf[,j] and cmbdf[i],
# and that pix are maintained.

# Check all generics like cbind and rbind

testthat::context("CMBDataFrame attributes")

a1 <- CMBDataFrame(nside = 1, ordering = "nested", coords = "spherical")
a2 <- CMBDataFrame(nside = 1, ordering = "nested", coords = "cartesian")
a3 <- CMBDataFrame(nside = 1, ordering = "nested")

b1 <- CMBDataFrame(nside = 1, ordering = "ring", coords = "spherical")
b2 <- CMBDataFrame(nside = 1, ordering = "ring", coords = "cartesian")
b3 <- CMBDataFrame(nside = 1, ordering = "ring")

testthat::test_that("coords attribute", {
  testthat::expect_equal(coords(a1), "spherical")
  testthat::expect_equal(coords(a2), "cartesian")
  testthat::expect_equal(coords(a3), NULL)
  testthat::expect_equal(coords(b1), "spherical")
  testthat::expect_equal(coords(b2), "cartesian")
  testthat::expect_equal(coords(b3), NULL)
})

testthat::test_that("ordering attribute", {
  testthat::expect_equal(ordering(a1), "nested")
  testthat::expect_equal(ordering(a2), "nested")
  testthat::expect_equal(ordering(a3), "nested")
  testthat::expect_equal(ordering(b1), "ring")
  testthat::expect_equal(ordering(b2), "ring")
  testthat::expect_equal(ordering(b3), "ring")
})

testthat::test_that("nside attribute", {
  testthat::expect_equal(nside(a1), 1)
  testthat::expect_equal(nside(a2), 1)
  testthat::expect_equal(nside(a3), 1)
  testthat::expect_equal(nside(b1), 1)
  testthat::expect_equal(nside(b2), 1)
  testthat::expect_equal(nside(b3), 1)
})

testthat::test_that("window attribute", {
  testthat::expect_equal(window(a1), NULL)
  testthat::expect_equal(window(a2), NULL)
  testthat::expect_equal(window(a3), NULL)
  testthat::expect_equal(window(b1), NULL)
  testthat::expect_equal(window(b2), NULL)
  testthat::expect_equal(window(b3), NULL)
})

testthat::test_that("is.CMBDataFrame TRUE", {
  testthat::expect_equal(is.CMBDataFrame(a1), TRUE)
  testthat::expect_equal(is.CMBDataFrame(a2), TRUE)
  testthat::expect_equal(is.CMBDataFrame(a3), TRUE)
  testthat::expect_equal(is.CMBDataFrame(b1), TRUE)
  testthat::expect_equal(is.CMBDataFrame(b2), TRUE)
  testthat::expect_equal(is.CMBDataFrame(b3), TRUE)
})



testthat::context("CMBDataFrame generics")

testthat::test_that("Square bracket operator", {
  testthat::expect_equal_to_reference(a1[1,], "references/square_bracket_a1[1,].rds")
  testthat::expect_equal_to_reference(a1[,1], "references/square_bracket_a1[,1].rds")
  testthat::expect_equal_to_reference(a1[1],  "references/square_bracket_a1[1].rds")
})


testthat::context("CMBDataFrame other arguments")

c1 <- CMBDataFrame(nside = 1, ordering = "nested", intensities = seq(1.01,1.12, by = 0.01))
c2 <- CMBDataFrame(nside = 1, ordering = "nested", win = list(CMBWindow(x = 0, y = 0, z = 1, r = 1),
                                                              CMBWindow(theta = c(0,1,1), phi = c(0,0,1))))

testthat::test_that("intensities argument", {
  testthat::expect_equal(c1$I, seq(1.01,1.12, by = 0.01))
})

testthat::test_that("window argument", {
  testthat::expect_equal(window(c2), list(CMBWindow(x = 0, y = 0, z = 1, r = 1),
                                          CMBWindow(theta = c(0,1,1), phi = c(0,0,1))))
})




a1 <- CMBDataFrame(nside = 1, ordering = "nested", coords = "spherical")
b1 <- cbind(a1, m = rep(1, 12))
c1 <- coords(a1, new.coords = "cartesian")
d1 <- coords(a1, new.coords = "spherical")
e1 <- coords(b1, new.coords = "cartesian")
f1 <- coords(b1, new.coords = "spherical")

a2 <- CMBDataFrame(nside = 1, ordering = "nested", coords = "cartesian")
b2 <- cbind(a2, m = rep(1, 12))
c2 <- coords(a2, new.coords = "cartesian")
d2 <- coords(a2, new.coords = "spherical")
e2 <- coords(b2, new.coords = "cartesian")
f2 <- coords(b2, new.coords = "spherical")

a3 <- CMBDataFrame(nside = 1, ordering = "nested")
b3 <- cbind(a3, m = rep(1, 12))
c3 <- coords(a3, new.coords = "cartesian")
d3 <- coords(a3, new.coords = "spherical")
e3 <- coords(b3, new.coords = "cartesian")
f3 <- coords(b3, new.coords = "spherical")

a4 <- CMBDataFrame(nside = 1, ordering = "nested")[,-1]
b4 <- cbind(a4, m = rep(1, 12))
c4 <- coords(a4, new.coords = "cartesian")
d4 <- coords(a4, new.coords = "spherical")
e4 <- coords(b4, new.coords = "cartesian")
f4 <- coords(b4, new.coords = "spherical")
