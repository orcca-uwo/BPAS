

# BPAS Polynomial System Solver

In this directory we supply a "portable" stand-alone polynomial system 
solver based on BPAS's regular chains.

* It is portable in the sense that it statically links to many 
  (but not all) required libraries. But, those non-static libraries
  should be available as shared libraries on most systems. 
  This binary statically links to the BPAS library, and so allows
  BPAS to be built, a solver built, and then BPAS built again with
  a different configuration (serial ON/OFF, generators ON/OFF). 

* It is stand-alone because it simply takes a file describing a polynomial
  system, solves it, and then returns the result. Such files should begin
  with a line of comma-separated variables describing the ambient space.
  Then, one line per polynomial in th einput system. For example, 
  see http://bpaslib.org/src/BPAS-Systems.tar.gz for a collection of
  many such files of polynomial systems. 

You may use the supplied Makefile to create the binary. However, note: 
* The BPAS library which is linked to this solver much be build with
  BPAS_WITH_MAPLE set to ON.

* You should set the LD_LIBRARY_PATH variable to include build: 
  ```export LD_LIBRARY_PATH=../../build/:$LD_LIBRARY_PATH``` 

* This solver requires an active Maple installation to communicate with
  at runtime. In particular, you will need to environment variables
  to point to your maple installation. In particular, LD_LIBRARY_PATH,
  CPLUS_INCLUDE_PATH, LIBRARY_PATH, and MAPLE. For example, if your Maple 
  installation is contained in /opt/maple2019 then you have the following:

```
export MAPLE=/opt/maple2019
export LIBRARY_PATH=$LIBRARY_PATH:$MAPLE/bin.X86_64_LINUX
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MAPLE/bin.X86_64_LINUX
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$MAPLE/extern/include
```

