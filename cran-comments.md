## Test environments

* local Windows 10 install, R 3.5.1
* local macOS "High Sierra" version 10.13.6, R 3.4.4
* Ubuntu 14.04 (on travis-ci) (devel, release and old release)
* winbuilder (devel and release)

## R CMD Check results
There were no ERRORS or WARNINGS.  
There were 2 NOTES:  


* Checking compiled code ... NOTE  
  File 'rcosmo/libs/x64/rcosmo.dll':  
    Found no calls to: 'R_registerRoutines',  
    'R_useDynamicSymbols'
    
    This note did not show up on winbuilder, Travis-CI, 
    or macOS. Followed all intructions: [here](https://stackoverflow.com/questions/42313373/r-cmd-check-note-found-no-calls-to-r-registerroutines-r-usedynamicsymbols) but could not remove this note.


* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Daniel Fryer <d.fryer@latrobe.edu.au>' 
    
    New submission    
    Possibly mis-spelled words in DESCRIPTION:    
    HEALPix (13:5, 14:15)    
    
    This word is not mis-spelled.


## Downstream dependencies
There are no downstream dependencies for this package

