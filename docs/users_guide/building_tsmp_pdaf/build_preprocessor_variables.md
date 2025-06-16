# TSMP-PDAF Preprocessor variables #

A number of preprocessor variables are set during the build process of
TSMP-PDAF in order to make the code behave in specific ways.

TSMP-PDAF environemnt variables are set in the machine-specific
TSMP-PDAF build-script:

```
bldsva/intf_DA/pdaf/arch/build_interface_pdaf.ksh
```

(precompiler:example)=
## Example for setting preprocessor variables ##

Before listing the preprocessor variables, we give an example for
setting them.

Preprocessor variables can be set in the TSMP-PDAF-build script
`TSMP/bldsva/intf_DA/pdaf/arch/build_interface_pdaf.ksh`.

A preprocessor variable `CPP_VAR` is set by adding a line of the form:
```bash
	cppdefs+=" ${pf}-DCPP_VAR "
```

If this line is added at the beginning of `configure_da()`, the
preprocessor variable is activated for all combinations of component
models.

Later in the code, preprocessor variable can be set only for specific combinations of component models. Here is an example source code from script `build_interface_pdaf.ksh`:
```bash
  if [[ $withCLM == "true" && $withCOS == "false" && $withPFL == "true" ]] ; then
     importFlags+=$importFlagsCLM
     importFlags+=$importFlagsOAS
     importFlags+=$importFlagsPFL
     importFlags+=$importFlagsDA
     cppdefs+=" ${pf}-DCOUP_OAS_PFL ${pf}-DMAXPATCH_PFT=1 "
     cppdefs+=" ${pf}-DOBS_ONLY_PARFLOW " # Remove for observations from both ParFlow + CLM
     if [[ $readCLM == "true" ]] ; then ; cppdefs+=" ${pf}-DREADCLM " ; fi
     if [[ $freeDrain == "true" ]] ; then ; cppdefs+=" ${pf}-DFREEDRAINAGE " ; fi
     libs+=$libsCLM
     libs+=$libsOAS
     libs+=$libsPFL
     obj+=' $(OBJCLM) $(OBJPF) '
  fi
```

The if-condition states that component models CLM and ParFlow are
used. This corresponds to the flag `-c clm-pfl-pdaf` in the
build-command.

The line `cppdefs+=" ${pf}-DOBS_ONLY_PARFLOW " # Remove for
observations from both ParFlow + CLM` sets
`OBS_ONLY_PARFLOW`. Analogously, any other preprocessor variable could
be set.


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

