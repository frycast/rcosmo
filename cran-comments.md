## Test environments

* local Windows 10 install, R 3.5.1
* Ubuntu 14.04 (on travis-ci) (devel, release and old release)
* winbuilder (devel and release)

## R CMD Check results
0 errors | 0 warnings | 2 notes

* Both NOTEs were as follows:  
  Checking compiled code ... NOTE  
  File 'rcosmo/libs/x64/rcosmo.dll':  
    Found no calls to: 'R_registerRoutines',  
    'R_useDynamicSymbols'

      Followed all intructions [here](https://stackoverflow.com/questions/42313373/r-cmd-check-note-found-no-calls-to-r-registerroutines-r-usedynamicsymbols) but could not remove this note.

## Downstream dependencies
There are no downstream dependencies for this package



