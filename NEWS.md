# rcosmo 1.0.0
This is the first release.


# rcosmo 1.01.0
* introduced depth_test to plot functions so that objects such
  as window boundaries and plot labels are not easily obscured
* more effective use of memory mapping to avoid reading data into
  memory unnecessarily, especially when using the win argument
  of CMBDataFrame function
* plot.CMBWindow now closes the polygon for neater visualisation
* fixed a small typo bug in nestSearch function
* updated links in rcosmo.R for Update and BugReports
* updated link in downloadCMBMap for data source with nside = 2048
* introduced check in nestSearch that log2 of nside must be integer
* fixed a small typo bug in qstat function
* new function geoAngle was added

