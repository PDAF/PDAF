# TSMP-PDAF Preprocessor variables #

A number of preprocessor variables are set during the build process of
TSMP-PDAF in order to make the code behave in specific ways.

TSMP-PDAF preprocessor variables are set in the TSMP2 CMake scripts
`BuildPDAF.cmake` and `BuildPDAFMODEL.cmake`.

(precompiler:example)=
## Example for setting preprocessor variables ##

Before listing the preprocessor variables, we give an example for
setting them.

Preprocessor variables can be set in the TSMP2 cmake scripts
`BuildPDAF.cmake` and `BuildPDAFMODEL.cmake`.

A preprocessor variable `CPP_VAR` is set by adding a line of the form:
```cmake
    list(APPEND PDAF_DEFS "-DCPP_VAR")
```

```{hint}
Check that the preprocessor variable is set for the right model
combinations. 
```

For example in 
```cmake
if(DEFINED OASIS_SRC)
  list(APPEND PDAF_DEFS "-Duse_comm_da")
  list(APPEND PDAF_DEFS "-DMAXPATCH_PFT=1")
  if(DEFINED PARFLOW_SRC)
    list(APPEND PDAF_DEFS "-DCOUP_OAS_PFL")
    list(APPEND PDAF_DEFS "-DOBS_ONLY_PARFLOW")
  endif()
  if(DEFINED eCLM_SRC)
    list(APPEND PDAF_DEFS "-DCLMFIVE")
  endif()
else()
  if(DEFINED CLM35_SRC)
    list(APPEND PDAF_DEFS "-DCLMSA")
  endif()
  if(DEFINED eCLM_SRC)
    list(APPEND PDAF_DEFS "-DCLMSA")
    list(APPEND PDAF_DEFS "-DCLMFIVE")
  endif()
endif()
```

`CLMSA` is defined iff 
- there is no coupling to ParFlow or ICON (`DEFINED OASIS_SRC` is
false)
- CLM is used in version 3.5 (`DEFINED CLM35_SRC` is true) or as eCLM
(`DEFINED eCLM_SRC` is true).

(precompiler:clmsa)=
## CLMSA ##

Preprocessor variable `CLMSA` is true if CLM-standalone is used in
TSMP-PDAF, i.e. no coupling to other component models (ParFlow,
atmospheric model).

`CLMSA` is used in many places in TSMP-PDAF, where CLM-standalone
specific code is introduced. This includes
- observation reading
- setting observation vector
- setting state vector
- communicator handling
- localized filters
- TSMP-PDAF-wrapper routines

(precompiler:parflow_stand_alone)=
## PARFLOW_STAND_ALONE ##

Preprocessor variable `PARFLOW_STAND_ALONE` is true if
ParFlow-standalone is used in TSMP-PDAF, i.e. no coupling to other
component models (CLM, atmospheric model).

It is used less frequently than [](precompiler:clmsa), only at code
places where the behavior of ParFlow-CLM-PDAF and ParFlow-PDAF should
differ.

(precompiler:obs_only_parflow)=
## OBS_ONLY_PARFLOW ##

Preprocessor variable `OBS_ONLY_PARFLOW` is true if observations in
TSMP-PDAF are of ParFlow-type.

This will remove unnecessary code during observation reading, when
ParFlow-CLM-PDAF is built, but no CLM-type observations are included.

(precompiler:obs_only_clm)=
## OBS_ONLY_CLM ##

Preprocessor variable `OBS_ONLY_CLM` is true if observations in
TSMP-PDAF are of CLM-type.

This will remove unnecessary code during observation reading, when
ParFlow-CLM-PDAF is built, but no ParFlow-type observations are
included.

## WATSAT3D ##

Preprocessor variable `WATSAT3D` is set in `common_build_interface.ksh`
in function `c_configure_clm`.

If it is turned on, the possibility of a read-in porosity is
implemented in CLM's `iniTimeConst.F90`.

Usage: Enforce consistent porosity in ParFlow and CLM.

## CLMFIVE ##

Currently only in feature branch `TSMP_pdaf-clm5`.

If defined, CLM5.0 is used.

If undefined, CLM3.5 is used.

This distinction is important in many parts of the TSMP-PDAF-wrapper
source code, when, f.e., function calls have changed from version 3 to
version 5 of CLM.

## OLD_TRUNCATE_SAT ##

If `OLD_TRUNCATE_SAT` is defined, the saturation truncation in
function `SaturationToPressure` (file
`problem_saturationtopressure.c`) is changed to an outdated default
(before October 2023).

Reference for the current default saturation truncation in Hung et al
(<https://doi.org/10.1029/2021WR031549>, Equation 11)

``` c++
	if(psdat[ips] <= (s_res + 0.003) ) psdat[ips] = s_res + 0.003;
```

`OLD_TRUNCATE_SAT` will change this to the outdated default:

``` c++
	if(psdat[ips] <= s_res) psdat[ips] = s_res + 0.01;
```

(precompiler:pdaf_debug)=
## PDAF_DEBUG ##

If `PDAF_DEBUG` is defined, PDAF's debugging output is turned on at
specific places in the code.

Currently implemented:
- `init_pdaf`
- `assimilate_pdaf`

Information on PDAF debugging:
<https://pdaf.awi.de/trac/wiki/PDAF_debugging>

