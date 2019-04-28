[![Travis-CI Build Status](https://travis-ci.org/frycast/rcosmo.svg?branch=master)](https://travis-ci.org/frycast/rcosmo) [![CRAN Versions](http://www.r-pkg.org/badges/version/rcosmo)](https://cran.r-project.org/web/packages/rcosmo/index.html) [![CRAN release dates](http://www.r-pkg.org/badges/version-ago/rcosmo)](https://cran.r-project.org/web/packages/rcosmo/index.html) [![CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/rcosmo)](https://cran.r-project.org/web/packages/rcosmo/index.html) [![CRAN downloads](http://cranlogs.r-pkg.org/badges/last-week/rcosmo)](https://cran.r-project.org/web/packages/rcosmo/index.html) [![Code coverage](https://codecov.io/gh/frycast/rcosmo/branch/master/graph/badge.svg)](https://codecov.io/github/frycast/rcosmo?branch=master)

# The `rcosmo` project

The Cosmic Microwave Background (CMB) is remnant electromagnetic radiation from the epoch of recombination. It is the most ancient important source of data about the early universe and the key to unlocking the mysteries of the Big Bang and the structure of time and space. Spurred on by a wealth of satellite data, intensive investigations in the past few years have resulted in many physical and mathematical results to characterize CMB radiation. An advanced R programming toolkit is needed to help statisticians perform CMB data analytics. 

`rcosmo` addresses various data processing and statistical analysis needs for the present generation of CMB experiments. These needs fall into the following broad categories:
+ Importing and transforming HEALPix data in convenient CMBDataFrames
+ Geometric tools
+ Statistical tools
+ Visualisation

The current version of `rcosmo` includes the following functionality:
+	Generation of a comprehensive data frame of CMB observations, which include HEALPix indices, metadata, CMB intensities and their
  corresponding spherical and/or cartesian coordinates, as well as polarization data
+	Window subsetting tools for investigating circular, convex and non-convex polygonal sub-regions on the sphere
+	Fast empirical covariance and variogram estimation
+	Implementation and analysis of spherical harmonics, spherical wavelets, etc
+	Various methods for CMB map visualization, such as interactive 3D full sky maps rendered with OpenGL, Mollweide projection and HEALPix boundary plotting
+ Spherical geometry tools such as shortest distance between two points, calculate spherical angles, shortest distance between a point and a region, etc



