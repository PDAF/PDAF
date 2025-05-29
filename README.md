
# PDAF (Parallel Data Assimilation Framework)

Copyright 2004-2025, Lars Nerger, Alfred Wegener Institute, Helmholtz Center
for Polar and Marine Research, Bremerhaven, Germany. 
For license information, please see the file LICENSE.txt.

For full documentation and tutorial, see: http://pdaf.awi.de 


## Note on PDAF V3.0

In the upgrade to PDAF V3.0 there are changes which make PDAF3
fully backward-compatible. If one has a code implemented for PDAF2,
one needs a few adaptions. For more information, see:
https://pdaf.awi.de/trac/wiki/PortingToPDAF3

## Introduction

PDAF is a unified framework for data assimilation with the aim
to reduce development effort for data assimilation while providing
high computational efficiency.

PDAF can be used
- to perform data assimilation with high-dimensional models,
- to assess data assimilation methods with small models,
- to teach ensemble data assimilation, and
- to develop new data assimilation methods and test them in a unified environment

PDAF provides
- A modular structure with separation-of-concerns for the model,
  the data assimilation methods, observation handling, and covariances
- A parallel infrastructure, using MPI and OpenMP, to implement a
  parallel data assimilation system based on existing numerical
  models (typically of, but not limited to, components of the Earth system). 
- A selection of state-of-the-art ensemble data assimilation algorithms
  based on the Kalman filter or nonlinear filters (see list below)
- A selection of 3D variational methods, both with parameterized
  and ensemble covariance matrix
- A structured approach to handle large sets of observation
  types
- Functionality for ensemble and observation diagnostics
- Functionality to generate synthetic observations for data
  assimilation studies (e.g. OSSEs)

The PDAF release provides also
- Tutorial codes demonstrating the implementation
- Code templates to assist in the implementation 
- Toy models fully implemented with PDAF for the study of data
  assimilation methods.
- Links to codes for existing model couplings with PDAF

## First Steps with PDAF

A good starting point for using PDAF is to run a tutorial example.
The web site  http://pdaf.awi.de/trac/wiki/FirstSteps 
provides hints on getting started with PDAF.

## PDAF tutorials
The page
https://pdaf.awi.de/trac/wiki/PdafTutorial
provides tutorial slide sets that explain the implementation
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
In order to build the PDAF library one needs a Fortran 2003 compiler,
'make'. PDAF is generally intended for parallel computing using MPI and
needs an MPI-library for compilation. In addition, the libraries BLAS
and LAPACK are required.

1. Choose a suitable include file for the make process and/or edit or
create one. See the directory make.arch/ for several provided include files.
For standard Linux systems using gfortran and OpenMPI using
linux_gfortran.h will most likely work.

2. Set the environment variable $PDAF_ARCH to the name of the include
file (without ending .h). Alternatively you can specify PDAF_ARCH on
the command line when running 'make' in step 3.

3. In the main directory execute 'make' at the prompt, e.g.
'make PDAF_ARCH=linux_gfortran'. This will
compile the sources and generate two library files that include the 
ensemble filter methods and the 3D-Var methods including solvers
in the directory lib/.  

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


## Test suite

The directory tests/ contains a set of implementations for consistency tests.
This is more for 'internal use'. We use these implementations to validate PDAF. 

You can run the scripts analogously to that in the directory tutorial/:
**runtests_offline.sh** or **runtests_online.sh** perform details tests
on options for the offline and online coupled implementations.
The script **test_runtime.sh** runs tests on a larger model grid. This is
usually used to assess the effect of code changes on the run time.


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
- **ENSRF** (Ensemble square-root filter  using serial observation processing and covariance localization)
       [J. Whitaker, T. Hamill, Mon. Wea. Rev. 130 (2002) 1913-1924]
- **EAKF** (Ensemble Adjustment Filter using serial observation processing and covariance localization)
       [J. Anderson, Mon. Wea. Rev. (2003) 634-642] 

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


PDAF is written in, mainly, Fortran 2003. 
The compilation and execution has been tested on the different systems ranging from
notebook computers to supercomputers, e.g.:
- Linux
- MacOS
- Cray CLE
- NEC Aurora
- Microsoft Windows 10 with Cygwin


## Contact Information

Please send comments, suggestions, or bug reports to pdaf@awi.de or use the Issues in Github.
