## Check area functions for CMBDataFrame and CMBWindow
cmbdf <- CMBDataFrame(nside = 1024, ordering = "nested",
                      coords = "cartesian")
disc <- CMBWindow(x = 1, y = 0, z = 0, r = 0.5)
m.disc <- CMBWindow(x = 1, y = 0, z = 0, r = 0.5, set.minus = TRUE)
polygon <- CMBWindow(phi = c(0, pi/4, pi/4, pi/5),
                     theta = c(pi/2, pi/2, pi/4, pi/2 - pi/20))
m.polygon <- CMBWindow(phi = c(0, pi/4, pi/4, pi/5),
                       theta = c(pi/2, pi/2, pi/4, pi/2 - pi/20),
                       set.minus = TRUE)
disc.cmb <- window(cmbdf, new.window = disc)
m.disc.cmb <- window(cmbdf, new.window = m.disc)
poly.cmb <- window(cmbdf, new.window = polygon)
m.poly.cmb <- window(cmbdf, new.window = m.polygon)

# Complements should add to 4*pi
isTRUE(all.equal(area(disc.cmb) + area(m.disc.cmb), 4*pi))
isTRUE(all.equal(area(poly.cmb) + area(m.poly.cmb), 4*pi))
isTRUE(all.equal(area(disc) + area(m.disc), 4*pi))
isTRUE(all.equal(area(polygon) + area(m.polygon), 4*pi))

# Area difference between CMBWindow and CMBDataFrame should
# be accounted for
isTRUE(all.equal(   area(poly.cmb) - area(polygon),
                  -(area(m.poly.cmb) - area(m.polygon)) ))
isTRUE(all.equal(   area(disc.cmb) - area(disc),
                  -(area(m.disc.cmb) - area(m.disc)) ))

# CMBWindow and CMBDataFrame should contain similar areas
# the tolerances here were chosen by inspection for this window
isTRUE(all.equal( area(poly.cmb), area(polygon), tolerance = 5e-3 ))
isTRUE(all.equal( area(m.poly.cmb), area(m.polygon), tolerance = 5e-5 ))
isTRUE(all.equal( area(disc.cmb), area(disc), tolerance = 1.5e-5 ))
isTRUE(all.equal( area(m.disc.cmb), area(m.disc), tolerance = 7e-7 ))
