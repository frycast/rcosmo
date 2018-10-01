## Test environments

* local Windows 10 install, R 3.5.1
* local macOS 10.13 "High Sierra" install, R 3.5.1 
* Ubuntu 14.04 (on travis-ci) (devel, release and old release)
* winbuilder (devel and release)

## R CMD Check results
There were no ERRORS or WARNINGS.  
There were 2 NOTES:  

* Checking dependencies in R code ... NOTE    
    Unexported object imported by a ':::' call: ‘geoR:::plot.variogram’

      This is necessary since we want to ensure that plot.variogram     
      is called without dispatch to another plot function,     
      but plot.variogram is not exported from geoR.    


* Checking compiled code ... NOTE  
  File 'rcosmo/libs/x64/rcosmo.dll':  
    Found no calls to: 'R_registerRoutines',  
    'R_useDynamicSymbols'
    
         Followed all intructions [here](https://stackoverflow.com/questions/42313373/r-cmd-check-note-found-no-calls-to-r-registerroutines-r-usedynamicsymbols) but could not remove this note.

    
* Found the following assignments to the global environment:    
  File ‘rcosmo/R/Statistics.R’:   
    assign(".temp.theta", NULL, pos = 1)    
    assign(".temp.theta", theta, pos = 1)    
    
    Fix this before submitting to CRAN


## Downstream dependencies
There are no downstream dependencies for this package


## TO DO BEFORE SUBMISSION

* rcosmo, Rcosmo, rCosmo are all placeholder docs
