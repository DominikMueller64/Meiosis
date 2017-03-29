## Test environments
* local ubuntu 14.04 install, R 3.3.2
* ubuntu 12.04 (on travis-ci), R 3.3.2
* ubuntu 12.04 (on travis-ci), R-release 
* ubuntu 12.04 (on travis-ci), R-devel 
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE for build with R-devel:

*checking compiled code ... NOTE

File ‘Meiosis/libs/Meiosis.so’:

  Found no calls to: ‘R_registerRoutines’, ‘R_useDynamicSymbols’

It is good practice to register native routines and to disable symbol search.

See ‘Writing portable packages’ in the ‘Writing R Extensions’ manual.


## Downstream dependencies
This is the first release.
