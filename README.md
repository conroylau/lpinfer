# linearprog

This is a WIP repository for the R modules for inference in linear programs.

##### install cplexAPI on Mac:
```
R CMD INSTALL --configure-args="--with-cplex-include=/Applications/CPLEX_Studio_Community129/cplex/include/ --with-cplex-lib=/Applications/CPLEX_Studio_Community129/cplex/lib/x86-64_osx/static_pic/" /Users/conroylam/Downloads/cplexAPI_1.3.6.tar.gz 
```
##### install Rcplex on Mac:
```
R CMD INSTALL --configure-args="PKG_CFLAGS='-m64 -fPIC' 
PKG_CPPFLAGS=-I/Applications/CPLEX_Studio_Community129/cplex/include 
PKG_LIBS='-L/Applications/CPLEX_Studio_Community129/cplex/lib/x86-64_osx/static_pic 
-lcplex -m64 -lm -lpthread -framework CoreFoundation -framework IOKit'" /Users/conroylam/Downloads/Rcplex
```
