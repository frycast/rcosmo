# rcosmo 1.1.3 (under development)

* Change to STRICT_R_HEADERS to maintain Rcpp compatibility.
* Fix typo in HPDataFrame documentation

# rcosmo 1.1.2

* No changes were made. The resubmission was due to archiving by CRAN.

# rcosmo 1.1.1

### Featured changes

 * Out-of-bounds array index bug fixed in minDist_internal1 and maxDist_internal1


# rcosmo 1.1.0

### Featured changes

 * New function geoAngle is added,
 * improves HPDataFrame inferral of assumedUniquePix 
   and healpixCentered attributes,
 * extends as.CMBDataFrame to work on HPDataFrame and 
   adds drop.coords parameter to as.CMBDataFrame,
 * more effective use of memory mapping to avoid reading data into
   memory unnecessarily, especially when using the win argument
   of CMBDataFrame function,
 * coords.HPDataFrame now behaves more like coords.CMBDataFrame,
   where leaving new.coords missing will return only the
   value of the attribute named coords,
 * HPDataFrame and nestSearch have the option to save dot products 
   for each observation with the nearest HEALPix pixel center,
 * Added function numeric2col to allow easy conversion of numeric
   vectors to colour schemes.
  
### Minor changes and bug fixes

 * Introduces depth_test to plot functions so that objects such
   as window boundaries and plot labels are not easily obscured,
 * plot.CMBWindow now closes the polygon for neater visualisation,
 * fixes a small typo bug in nestSearch function,
 * updated links in rcosmo.R for Update and BugReports,
 * updated link in downloadCMBMap for data source with nside = 2048,
 * introduces check in nestSearch that log2 of nside must be integer,
 * fixes a small typo bug in qstat function,
 * bug fix HPDataFrame producing unique pixel indices when,
   delete.duplicates = TRUE,
 * coords(x) <- NULL now sets coordinates to NULL rather than
   producing an error,
 * introduces healpixCentered attribute getter function,
 * introduces option to save row indices of duplicated 
   pixel indices for HPDataFrame when delete.duplicates
   is set to TRUE.



# rcosmo 1.0.0
This is the first release.



   
   


