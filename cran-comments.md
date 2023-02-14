## Test environments
* local ubuntu 20.04 install, R 4.2.1
* win-builder (devel, release and old-release)
* r-hub (windows-server, Fedora Linux R-devel, clang, gfortran, solaris)

## R CMD check results

There were no ERRORs or WARNINGs. 
1 NOTE:
- there are ::: calls to the package's namespace in its code. I could not fix this note as it is used as argument in `substitute` function, and I do not want to export the function.