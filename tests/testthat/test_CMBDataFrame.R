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


# -------------------------------------------------
testthat::context("CMBDataFrame attributes")

a1 <- CMBDataFrame(nside = 2, ordering = "nested", coords = "spherical")
a2 <- CMBDataFrame(nside = 2, ordering = "nested", coords = "cartesian")
a3 <- CMBDataFrame(nside = 2, ordering = "nested")

b1 <- CMBDataFrame(nside = 2, ordering = "ring", coords = "spherical")
b2 <- CMBDataFrame(nside = 2, ordering = "ring", coords = "cartesian")
b3 <- CMBDataFrame(nside = 2, ordering = "ring")

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
  testthat::expect_equal(nside(a1), 2)
  testthat::expect_equal(nside(a2), 2)
  testthat::expect_equal(nside(a3), 2)
  testthat::expect_equal(nside(b1), 2)
  testthat::expect_equal(nside(b2), 2)
  testthat::expect_equal(nside(b3), 2)
})

testthat::test_that("window attribute", {
  testthat::expect_equal(window(a1), NULL)
  testthat::expect_equal(window(a2), NULL)
  testthat::expect_equal(window(a3), NULL)
  testthat::expect_equal(window(b1), NULL)
  testthat::expect_equal(window(b2), NULL)
  testthat::expect_equal(window(b3), NULL)
})

testthat::test_that("pix attribute", {
  testthat::expect_equal(pix(a1), seq(1,48))
  testthat::expect_equal(pix(a2), seq(1,48))
  testthat::expect_equal(pix(a3), seq(1,48))
  testthat::expect_equal(pix(b1), seq(1,48))
  testthat::expect_equal(pix(b2), seq(1,48))
  testthat::expect_equal(pix(b3), seq(1,48))
})

testthat::test_that("is.CMBDataFrame TRUE", {
  testthat::expect_equal(is.CMBDataFrame(a1), TRUE)
  testthat::expect_equal(is.CMBDataFrame(a2), TRUE)
  testthat::expect_equal(is.CMBDataFrame(a3), TRUE)
  testthat::expect_equal(is.CMBDataFrame(b1), TRUE)
  testthat::expect_equal(is.CMBDataFrame(b2), TRUE)
  testthat::expect_equal(is.CMBDataFrame(b3), TRUE)
})





# -------------------------------------------------
testthat::context("CMBDataFrame other arguments")

c1 <- CMBDataFrame(nside = 2, ordering = "nested", I = seq(1.01,1.48, by = 0.01))
c2 <- CMBDataFrame(nside = 2, ordering = "nested", win = list(CMBWindow(x = 0, y = 0, z = 1, r = 1), CMBWindow(theta = c(0,1,1), phi = c(0,0,1))))
c3 <- CMBDataFrame(nside = 2, ordering = "nested", I = seq(1.01,1.48, by = 0.01), spix = c(1,5,9,7))
c4 <- CMBDataFrame(nside = 2, ordering = "nested", I = seq(1.01,1.48,by=0.01), sample.size = 48)
c5 <- CMBDataFrame(nside = 2, ordering = "nested", I = seq(101,148), sample.size = 5)

testthat::test_that("intensities argument", {
  testthat::expect_equal(c1$I, seq(1.01,1.48, by = 0.01))
})

testthat::test_that("window argument", {
  testthat::expect_equal(window(c2), list(CMBWindow(x = 0, y = 0, z = 1, r = 1),
                                          CMBWindow(theta = c(0,1,1), phi = c(0,0,1))))
})

testthat::test_that("spix argument", {
  testthat::expect_equal(pix(c3), c(1,5,7,9))
  testthat::expect_equal(c3, c1[c(1,5,7,9),])
})

testthat::test_that("sample.size argument", {
  testthat::expect_equal(pix(c4), 1:48)
  testthat::expect_equal(c4$I, seq(1.01,1.48,by=0.01))
  testthat::expect_equal(all(c5$I - 100 == pix(c5)), TRUE)
})





# -------------------------------------------------
testthat::context("Passing CMBDF to CMBDF")

d1 <- CMBDataFrame(nside = 2, ordering = "nested", spix = c(2,3,4,11))
d2 <- CMBDataFrame(nside = 2, ordering = "ring", spix = c(2,3,4,11))
d3 <- CMBDataFrame(nside = 2, ordering = "nested", I = seq(1.01,1.48, by = 0.01), spix = c(2,3,4,11))

testthat::test_that("change coords", {
  testthat::expect_equal(a2, CMBDataFrame(a1, coords = "cartesian"))
  testthat::expect_equal(a1, CMBDataFrame(a2, coords = "spherical"))
  testthat::expect_equal(a1, CMBDataFrame(a3, coords = "spherical"))
  testthat::expect_equal(a2, CMBDataFrame(a3, coords = "cartesian"))
  testthat::expect_equal(b2, CMBDataFrame(b1, coords = "cartesian"))
  testthat::expect_equal(b1, CMBDataFrame(b2, coords = "spherical"))
  testthat::expect_equal(b1, CMBDataFrame(b3, coords = "spherical"))
  testthat::expect_equal(b2, CMBDataFrame(b3, coords = "cartesian"))
})

testthat::test_that("change coords on subset CMBDF", {
  testthat::expect_equal(a2[c(2,3,4,11),], CMBDataFrame(d1, coords = "cartesian"))
  testthat::expect_equal(a1[c(2,3,4,11),], CMBDataFrame(d1, coords = "spherical"))
  testthat::expect_equal(b2[c(2,3,4,11),], CMBDataFrame(d2, coords = "cartesian"))
  testthat::expect_equal(b1[c(2,3,4,11),], CMBDataFrame(d2, coords = "spherical"))
  testthat::expect_equal(a1[c(6,8,10,48),], CMBDataFrame(a2[c(6,8,10,48),], coords = "spherical"))
  testthat::expect_equal(a2[c(6,8,10,48),], CMBDataFrame(a1[c(6,8,10,48),], coords = "cartesian"))
})


testthat::test_that("use spix on subset CMBDF", {
  testthat::expect_equal(c1[c(3,11),], CMBDataFrame(d3, spix = c(3,11)))
})




# -------------------------------------------------
testthat::context("Reading to CMBDF from FITS")








# -------------------------------------------------
testthat::context("CMBDataFrame generics")

testthat::test_that("Square bracket operator", {
  testthat::expect_equal_to_reference(a1[1,], "references/square_bracket_a1[1,].rds")
  testthat::expect_equal_to_reference(a1[,1], "references/square_bracket_a1[,1].rds")
  testthat::expect_equal_to_reference(a1[1],  "references/square_bracket_a1[1].rds")
})




# -------------------------------------------------
testthat::context("CMBDataFrame helpers")

sp.nest <- c(1,9,32,45)
sp.ring <- nest2ring(nside = 2, pix = sp.nest)
e1 <- CMBDataFrame(nside = 2, ordering = "nested", coords = "cartesian", spix = sp.nest)
e2 <- CMBDataFrame(nside = 2, ordering = "ring", coords = "cartesian", spix = sp.ring)


testthat::test_that("Change of ordering scheme", {
  testthat::expect_equal(ordering(a1, new.ordering = "ring"), b1)
  testthat::expect_equal(ordering(a2, new.ordering = "ring"), b2)
  testthat::expect_equal(ordering(a3, new.ordering = "ring"), b3)
  testthat::expect_equal(ordering(b1, new.ordering = "nested"), a1)
  testthat::expect_equal(ordering(b2, new.ordering = "nested"), a2)
  testthat::expect_equal(ordering(b3, new.ordering = "nested"), a3)
  testthat::expect_equivalent(a2[sp.nest,], b2[sp.ring,])
  testthat::expect_equal(ordering(e1, new.ordering = "ring"), e2)
  testthat::expect_equal(ordering(e2, new.ordering = "nested"), e1)
})
