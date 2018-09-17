## Test environments

* local Windows 10 install, R 3.5.1
* <insert e.g., Ubuntu 12.04 (on travis-ci), R 3.5.1>
* <insert e.g., winbuilder (devel and release)>

## R CMD Check results
0 errors | 0 warnings | 2 notes

* Both NOTEs were as follows:  
  Checking compiled code ... NOTE  
  File 'rcosmo/libs/x64/rcosmo.dll':  
    Found no calls to: 'R_registerRoutines',  
    'R_useDynamicSymbols'

