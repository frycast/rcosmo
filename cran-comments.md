## Update
This is an update from version 1.0.0 to 1.01.0. See NEWS.md for 
the list of changes.


## Test environments

* local Windows 10 install, R 3.5.3
* Ubuntu 16.04.6 LTS (on travis-ci) (R devel, 3.6.0 and 3.5.3)
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


## Other comments
The manuscript describing the methods in this package is under preparation

