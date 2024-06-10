# TSMP-PDAF: Pre-patched PDAF

This fork of PDAF contains the PDAF changes ("patches") for using it
with TSMP (both TSMP1 <https://github.com/HPSCTerrSys/TSMP> and TSMP2
<https://github.com/HPSCTerrSys/TSMP2>).

1. the PDAF-library with slight changes in
   1. `make.arch/`: changed include files, in particular added the
include file `cmake.h` for TSMP2-PDAF
   2. `src/`: debug output of the observation ensemble for EnKF/LEnKF
2. `interface/` the interface routines for TSMP-PDAF

*Explanation of the term "pre-patched"*: For TSMP1-PDAF, the source
code of PDAF is changed during the build-process. F.e. PDAF source
code files may be replaced by corresponding files from the
TSMP1-repository.  The process of changing the source code of PDAF is
called "patching". For TSMP2-PDAF, there is no "patching" of
PDAF. Instead, TSMP2-PDAF clones the source code of this fork of
PDAF. This fork already includes all needed TSMP-related changes,
i.e. it is "pre-patched".

## Updating the pre-patched PDAF

The pre-patched PDAF is updated with three types of changes.

1. Developments from TSMP1-PDAF.
   - *Remark*: As long as TSMP1 is actively developed, PDAF-related
     development is carried out in TSMP1 and then copied to this
     repository. Once TSMP1 will be deprecated, the TSMP2-PDAF
     interface will be developed in this repository.
   - *Technical explanation*: Changes in the TSMP-PDAF interface are
     implemented in the repo `TSMP` and then copied (and re-commited)
     to `tsmp-pdaf-patched`, the default branch of this
     repository. These changes become part of TSMP2, when the
     pre-patched PDAF is cloned as described in TSMP2's
     README. Version tags (<https://github.com/HPSCTerrSys/pdaf/tags>)
     ensure stability and developments become part of the default
     TSMP2-build once the version tags are updated in the README.

2. Version updates from PDAF
   - *Technical explanation*: PDAF version updates are pulled from
     PDAF's branch `master` into the branch `tsmp-pdaf-patched`. For a
     PDAF version, f.e. `PDAF_V2.2.1` in the `master`, a corresponding
     version tag is defined, f.e. `PDAF_V2.2.1-tsmp` in
     `tsmp-pdaf-patched`.

3. TSMP-related changes in the PDAF library
   -  *Remark*: TSMP-related changes in the PDAF library should be
   avoided. If implemented, these changes should be described in this
   README.
   - *Technical explanation*: TSMP-related changes in the PDAF library
   are implemented directly in `tsmp-pdaf-patched`.


# PDAF (Parallel Data Assimilation Framework)

Copyright 2004-2024, Lars Nerger, Alfred Wegener Institute, Helmholtz Center
for Polar and Marine Research, Bremerhaven, Germany. 
For license information, please see the file LICENSE.txt.

For full documentation and tutorial, see: http://pdaf.awi.de 


## Introduction

PDAF is a framework for data assimilation.
PDAF can be used to assess data assimilation methods with small models,
to perform real data assimilation with high-dimensional models, and to
teach ensemble data assimilation. 

PDAF provides 
- A parallel infrastructure, using MPI and OpenMP, to implement a
  parallel data assimilation system based on existing numerical
  models (typically of components of the Earth system). 
- A selection of ensemble data assimilation algorithms based on
  the Kalman filter or nonlinear filters (see list below)
- A selection of 3D variational methods, both with parameterized
  and ensemble covariance matrix
- Functions for ensemble diagnostics
- Functionality to generate synthetic observations for data
  assimilation studies (e.g. OSSEs)

The PDAF release provides also
- Tutorial codes demontrating the implementation
- Code templates for assist in the implementation
- Toy models fully implemented with PDAF for the study of data
  assimilation methods.
- Model bindings for using PDAF with different models 


## First Steps with PDAF

A good starting point for using PDAF is to run a tutorial example.
The directory /tutorial contains files demonstrating the application 
of PDAF with a simple 2-dimensional example.

The web site  http://pdaf.awi.de/trac/wiki/FirstSteps 
provides hints on getting started with PDAF and 
  https://pdaf.awi.de/trac/wiki/PdafTutorial
holds the tutorial slide sets that explain the implementation
steps and how to compile and run the tutorial examples. 


## Models

The directory models/ contains toy models that are fully implemented
with PDAF. These models can be used to assess the behavior of different
assimilation algorithms. 
- **lorenz96/**
  This directory contains the Lorenz-96 model, which is a widely used
  model to assess data assimilation methods. Provided is a full
  implementation of the data assimilation with various filters and options.
  This model can be configured to have a sufficiently large state
  dimension to test low-rank filter algorithms.
- **lorenz63/**
  This directory contains the Lorenz-63 model, which is a classical
  3-variable model with chaotic dynamics. Provided is a full
  implementation of the data assimilation with various filters and options.
  The small state dimension and nonlinear dynamics make it a suitable
  test case for the standard particle filter (PF).
- **lorenz05b/**
  This directory contains the model Lorenz-2005 model II. Provided is a full
  implementation of the data assimilation with various filters and options.
- **lorenz05c/**
  This directory contains the two-scale model Lorenz-2005 model III.
  Provided is a full implementation of the data assimilation with various
  filters and options.

Instructions on how to run data assimilation experiments with these models
are provided on the PDAF web site.


## Installation of the PDAF library

The PDAF library will be automatically built when compiling a tutorial case
or one of the models. However, one can also separately build the library.
In order to build the PDAF library you need a Fortran 90 compiler, and
'make'

1. Choose a suitable include file for the make process and/or edit
one. See the directory make.arch/ for several provided include files.
There are include files for compilation with and without MPI.

Note: PDAF is generally intended for parallel computing using MPI.
However, it can be compiled for serial computing. To compile PDAF
for this case, a simplified MPI header file is included und should be
in the include path. In addition, a dummy implementation of MPI, which
behaves like MPI in the single-process case, is provided in the
directory nullmpi/. For the serial case, this file should also be
compiled and linked when PDAF is linked to a program.

2. Set the environment variable $PDAF_ARCH to the name of the include
file (without ending .h). Alternatively you can specify PDAF_ARCH on
the command line when running 'make' in step 3.

3. cd into the directory src/ and execute 'make' at the prompt. This will
compile the sources and generate a library file that includes the 
ensemble filter methods in the directory lib/. 
To generate the PDAF library including the 3D-Var methods and the 
solvers from the external libraries in /external/ execute
'make pdaf-var' at the prompt.
 


## Test suite

The directory testsuite/ contains another set of example implementations.
This is more for 'internal use'. We use these implementations to validate PDAF. 
The model is trivial: At each time step simply the time step size is added 
to the state vector. In this example all available filters are implemented.


## Verifying your installation 

The tutorial implementations can be verified as follows:

You can run the script
**runtests.sh** in the main tutorial directory **tutorial**.

This script will compile and run all tutorial implementations. Afterwards
the outputs at the final time step are checked against reference outputs
from the directory verification using a Python script. You can also compare
 the output files like out.online_2D_parallelmodel with reference files.
(Note: The script runtests.sh uses the generic compile definitions for
Linux with the gfotran compiler. For other systems, you might need to
change the settings for the make definitions files).


The testsuite also provides a functionality for verification:

Using 'make' one can run test cases for the verification which are
compared to reference outputs provided in the sub-directories
of the directory  testsuite/tests_dummy1D for different computers
and compilers. In particular the online case dummymodel_1D and the
offline test offline_1D can be run. Scripts for serial (non-parallel)
execution as well as example scripts for running parallel test jobs on
computers with SLURM or PBS batch systems are provided.

An installation of PDAF can be verified using the test suite as follows:
1. prepare the include file in make.arch
2. cd to testsuite/src
3. Build and execute the online experiments:
   'make pdaf_dummy_online' and
   'make test_pdaf_online > out.test_pdaf_online'
4. Build and execute the offline experiments:
   'make pdaf_dummy_online' and
   'make test_pdaf_offline > out.test_pdaf_offline'
6. Check the files out.test_pdaf_online and out.test_pdaf_offline
   At the end of the file, you see a list of Checks done using
   a Python script. Here the outputs are compared with reference 
   outputs produced with gfortran and MacOS.
   You can also diff the files to corresponding files in one of the
   example-directories in ../tests_dummy1D. Here, also reference
   output files, like output_lestkf0.dat are stored.


## Data Assimilation Algorithms 

The filter algorithms in PDAF are:

**Filters with global analysis step**
- **EnKF** (The classical perturbed-observations Ensemble Kalman filter)
       [G. Evensen, J. Geophys. Res. 99 C5 (1994) 10143-10162,
        G. Burgers et al., Mon. Wea. Rev. 126 (1998) 1719-1724]
- **ESTKF** (Error Subspace Transform Kalman filter)
       [L. Nerger et al. Mon. Wea. Rev. 140 (2012) 2335-2345, doi:10.1175/MWR-D-11-00102.1]
- **ETKF** (Ensemble Transform Kalman filter)
       [C. H. Bishop et al. Mon. Wea. Rev. 129 (2001) 420-436]
       The implementation in PDAF follows that described for the LETKF, but as a global filter. 
- **SEIK** (Singular "Evolutive" Interpolated Kalman) filter
       This is the full ensemble variant of the SEIK
       (Singular "Interpolated" Extended Kalman) filter.
       [SEIK: D.-T. Pham et al., C. R. Acad Sci., Ser. III, 326 (1009)
        255-260, for the SEIK variant in PDAF see L. Nerger et al.,
        Tellus 57A (2005) 715-735, doi:10.3402/tellusa.v57i5.14732]	
- **SEEK** (Singular "Evolutive" Extended Kalman) filter
       [D.-T. Pham et al., J. Mar. Syst. 16 (1998) 323-340] 
- **NETF** (Nonlinear Ensemble Transform Filter)
       [J. Toedter, B. Ahrens, Mon. Wea. Rev. 143 (2015) 1347-1367, doi:10.1175/MWR-D-14-00108.1]
- **PF** (Particle filter with importance resampling)
       [see S. Vetra-Carvalho et al., Tellus A 70 (2018) 1445364, doi:10.1080/16000870.2018.1445364]

**Filters with localized analysis step**
- **LESTKF** (Local Error Subspace Transform Kalman filter)
       [L. Nerger et al. Mon. Wea. Rev. 140 (2012) 2335-2345, doi:10.1175/MWR-D-11-00102.1]
- **LETKF** (Local Ensemble Transform Kalman filter)
       [B. R. Hunt et al., Physica D 230 (2007) 112-126, doi:10.1016/j.physd.2006.11.008]
- **LSEIK** (Local Singular "Evolutive" Interpolated Kalman) filter
       [L. Nerger et al., Oce. Dyn. 56 (2006) 634-649, doi:10.1007/s10236-006-0083-0]
- **LEnKF** (The classical perturbed-observations Ensemble Kalman filter with localization)
       [G. Evensen, J. Geophys. Res. 99 C5 (1994) 10143-10162,
        G. Burgers et al., Mon. Wea. Rev. 126 (1998) 1719-1724]
- **LNETF** (Nonlinear Ensemble Transform Filter with observation localization)
       [J. Toedter, B. Ahrens, Mon. Wea. Rev. 143 (2015) 1347-1367, doi:10.1175/MWR-D-14-00108.1]
- **LKNETF** (Local Kalman-nonlinear ensemble transform filter)
       [L. Nerger, Q. J. R. Meteorol Soc., 148 (2022) 620-640, doi:10.1002/qj.4221]

All filter algorithms are fully parallelized with MPI and optimized. The local filters 
(LSEIK, LETKF, LESTKF, LNETF, LKNETF) are in addition parallelized using OpenMP.

**Smoother extensions** are included for the filters ESTKF/LESTKF, ETKF/LETKF, EnKF, NETF/LNETF.

**3D-Var methods**

Next to the ensemble filter methods, different 3D-Var methods are provided:
- **parameterized 3D-Var**
- **3D ensemble** Var with ensemble transformation by ESTKF or LESTKF
- **hybrid 3D-Var** with ensemble transformation by ESTKF or LESTKF

The 3D-Var methods are implemented as incremental 3D-Var schemes following
Bannister, Q. J. Royal Meteorol. Soc., 143 (2017) 607-633, doi:10.1002/qj.2982. 


PDAF is written in Fortran (mainly Fortran 90 with some features from Fortran 2003). 
The compilation and execution has been tested on the different systems ranging from
notebook computers to supercomputers, e.g.:
- Linux
- MacOS
- Cray CLE
- NEC Super-UX
- Microsoft Windows 10 with Cygwin


## Contact Information

Please send comments, suggestions, or bug reports to pdaf@awi.de
