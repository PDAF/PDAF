(enkfpf)=
# Control file `enkfpf.par` #

An additional file named `enkfpf.par` needs to be present in the
TSMP-PDAF run directory.

`enkfpf.par` is read by the routine `read_enkfpar` in
`model/common/read_enkfpar.c`.

`enkfpf.par` contains information about the different model components
(ParFlow, CLM and data assimilation) like, e.g., the timing
information of the models, the prefixes for the ParFlow/ CLM input
files, etc.

This input is structured in four sections: 
- [`[PF]`](enkfpf:pf) which holds information about ParFlow, 
- [`[CLM]`](enkfpf:clm) which holds information about CLM 
- [`[COSMO]`](enkfpf:cosmo) which holds information about COSMO
- [`[DA]`](enkfpf:da) which holds information about the data assimilation process

An example for `enkfpf.par` is given below.  Note that the sequence of
the entries within these four categories can be changed if desired.

``` text
[PF]
problemname        = ""
nprocs             =
starttime          =
dt                 =
endtime            =
simtime            =
updateflag         =
gwmasking          =
paramupdate        =
paramupdate_frequency =
dampingfactor_param =
dampingfactor_state =
damping_switch_sm =
aniso_perm_y       =
aniso_perm_z       =
aniso_use_parflow  =
printensemble      =
t_printensemble    =
printstat          =
paramprintensemble =
paramprintstat     =
olfmasking         =
olfmasking_param   =
olfmasking_depth   =

[CLM]
problemname = ""
nprocs      =
update_swc  =
update_texture  =
update_T  =
print_swc   =
print_et   =
statevec_allcol =
statevec_colmean =
statevec_only_active =
statevec_max_layer =
t_printensemble =
watmin_switch =

[COSMO]
nprocs      =
dtmult      =

[DA]
outdir = ""
nreal =
startreal =
total_steps =
tstartcycle =
da_interval       =
da_interval_final =
flexible_da_interval =
stat_dumpoffset   =
point_obs =
crns_flag =
da_crns_depth_tol =
print_obs_index =
```

In the following the individual entries of `enkfpf.par` are described:
(enkfpf:pf)=
## [PF] ##

### PF:problemname ###

`PF:problemname`: (string) Problem prefix for ParFlow.

(enkfpf:pf:nprocs)=
### PF:nprocs ###

`PF:nprocs`: (integer) Number of processors per ParFlow instance. Must
match with the specifications in the `*.pfidb` input.

This number of processors specifies a subset of the processors
available in a single `COMM_model`. Each `COMM_model` contains
`npes_world / nreal` processes.

### PF:starttime ###

`PF:starttime`: (real) ParFlow start time. Must match with the
specifications in the `*.pfidb` input (`TimingInfo.StartTime`).

### PF:dt ###

`PF:dt`: (real) Inverse factor multiplied to CLM and COSMO time
steps. For coupled simulations, the quotient of CLM/COSMO time steps
and this factor should be in agreement with the ParFlow time step
(`TimingInfo.BaseUnit` in `*.pfidb` input files).

For CLM-ParFlow, `PF:dt` should be set as
`PF:dt = (dtime / TimingInfo.BaseUnit)`, where `dtime` is
the time delta specified in the CLM namelist file (`lnd.stdin`).

For COSMO-ParFlow, `PF:dt` should be set as
`PF:dt = (dt_cos / TimingInfo.BaseUnit) * COSMO:dtmult`.
TODO: Document COSMO time step setting `dt_cos`.

It is implicitly assumed that ParFlow and CLM calculate the same
number of time steps between two PDAF-calls (this number of time steps
is specified in `DA:da_interval`).

For CLM-standalone simulations `PF:dt` determines the time unit of
`DA:da_interval` for CLM simulations as `(dtime / PF:dt)`. For
CLM-standalone simulations, matching with the ParFlow time step is not
an issue and it makes sense to choose `PF:dt==1.0`, such that
`DA:da_interval` will specify the number of CLM time steps.

#### Examples for PF:dt ####

**Example 1**, CLMSA: `PF:dt==1` means that CLM input `dtime` is the unit of
`DA:da_interval`. For the default `dtime` of half an hour (1800
seconds), `DA:da_interval==48` would specify daily observations.

**Example 2**, CLMSA: For `dtime==1800` and `PF:dt==0.5`, the unit of
`DA:da_interval` is 1 hour - the standard time unit of
ParFlow. `DA:da_interval==24` would specify daily observations.

**Example 3**, FallSchoolCase, CLM-ParFlow-PDAF: The FallSchoolCase
has the following inputs concerning time stepping
- CLM: `dtime=3600`
- ParFlow: `pfset TimingInfo.BaseUnit 1.0`
- `enkfpf.par`: (1) `PF:dt==1.0`, (2) `DA:da_interval==1`

Both, CLM and ParFlow use a time step of 1 hour. Therefore,
`PF:dt==1.0` is the correct choice. In this setup, `DA:da_interval==1`
specifies hourly observations.

**Example 4**: eCLM-PDAF, open loop: For an eCLM open loop runs it makes
sense to choose `DA:da_interval` equal to `PF:simtime`, such that only
one "assimilation" cycle is computed. Additionally, it makes sense to
choose `PF:dt==1`, which means that `DA:da_interval` specifies
`tsclm`, the number of time steps that eCLM will compute.

\begin{align*}
\mathtt{tsclm} &= \frac{\mathtt{DA:da\_interval}}{\mathtt{PF:dt}}
\end{align*}


### PF:endtime (deprecated) ###

Deprecated. Sets `PF:simtime`.

### PF:simtime ###

`PF:simtime`: (real) Total simulation time (in terms of ParFlow
timing). 

For all simulations (including CLMSA): `PF:simtime` is used in
conjunction with `DA:da_interval` to determine `total_steps`, the
total number of iterations of the main data assimilation loop. Each
iteration of the main data assimilation loop consists of one forward
simulation and one data assimilation step.

\begin{align*}
\mathtt{total\_steps} &= \frac{\mathtt{PF:simtime}}{\mathtt{DA:da\_interval}}
\end{align*}

For ParFlow simulations, `PF:simtime` must match with the
specifications in the `*.pfidb` input: `PF:simtime` must correspond to
`TimingInfo.StopTime` MINUS `TimingInfo.StartTime`!


### PF:updateflag ###

`PF:updateflag`: (integer) Type of state vector update in ParFlow.

-   1: Assimilation of pressure data. State vector consists of
    pressure values (groundwater masking or a mixed state vector is
    introduced through `PF:gwmasking`) and is directly updated with
    pressure observations.

-   2: Assimilation of soil moisture data. State vector consists of
    soil moisture content values and is updated with soil moisture
    observations. The updated soil moisture content is transformed
    back to pressure via the inverse van Genuchten functions.

-   3: Assimilation of soil moisture data. State vector consists of
    soil moisture content and pressure values. Soil moisture data are
    used to update pressure indirectly.

### PF:gwmasking ###

`PF:gwmasking`: (integer) Groundwater masking for assimilation of
pressure data (`updateflag=1`) in ParFlow.

-   0: No groundwater masking.

-   1: Groundwater masking using saturated cells only.

-   2: Groundwater masking using mixed state vector.

For assimilating SM data (`updateflag=2`), the following masking
options are included:

-   0: No groundwater masking.

-   1: Groundwater masking using **unsaturated** cells only.

### PF:paramupdate ###

`PF:paramupdate`: (integer) Flag for parameter update

-   0: No parameter update

-   1: Update of saturated hydraulic conductivity

-   2: Update of Mannings coefficient

- 3: Update porosity

- 4: Update van Genuchten parameters

- 5: Update hydraulic conductivity and porosity

- 6: Update hydraulic conductivity and van Genuchten parameters

- 7: Update porosity and van Genuchten parameters

- 8 Update hydraulic conductivity, porosity and van Genuchten parameters


### PF:paramupdate_frequency ###
`PF:paramupdate_frequency`: (integer) Frequency of parameter
updates. Default: `1`

For each assimilation cycle it is checked whether `tcycle mod
paramupdate_frequency == 0`. When this happens, a parameter update is
applied. For the default of `1`, each assimilation cycle contains a
parameter update.

(enkfpf:pf:dampingfactor_param)=
### PF:dampingfactor_param ###
`PF:dampingfactor_param`: (real) Damping factor for parameter
updates. Default `1.0`.

The damping factor should be chosen between `0.0` and `1.0`, where the
inputs yield the following behavior:
- `0.0`: No parameter update.
- `1.0`: Parameter update without damping (default)

General equation for the damped, updated parameter vector

\begin{align*}
p_{update, damped} &= p + \mathtt{dampingfactor\_param} \cdot ( p_{update} - p )
\end{align*}

where $p$ is the parameter vector that is part of the state vector
before assimilation and $p_{update}$ is $p$ after the assimilation.

(enkfpf:pf:dampingfactor_state)=
### PF:dampingfactor_state ###
`PF:dampingfactor_state`: (real) Damping factor for state
updates. Default `1.0`.

Remark: Currently implemented for `updateflag==1`.

The damping factor should be chosen between `0.0` and `1.0`, where the
inputs yield the following behavior:
- `0.0`: No state update.
- `1.0`: State update without damping (default)

General equation for the damped, updated state vector

\begin{align*}
x_{update, damped} &= x + \mathtt{dampingfactor\_state} \cdot ( x_{update} - x )
\end{align*}

where $x$ is the state vector (without parameters) before assimilation
and $x_{update}$ is the state vector after the assimilation.

### PF:damping_switch_sm ###
`PF:damping_switch_sm`: (integer) Switch for applying damping factor
for state updates to soil moisture. Default `1` (damping factor
applied). Only applies for `PF:gwmasking==2`.

- `0`: State damping not applied to soil moisture
- `1`: State damping applied to soil moisture

General state damping is turned on by setting `PF:dampingfactor_state`
to a positive value smaller than `1.0`.

### PF:aniso_perm_y ###

`PF:aniso_perm_y`: (real) Anisotropy factor of saturated hydraulic
conductivity in y-direction. Only used when hydraulic conductivity is
updated (`[PF]paramupdate = 1`). Important: Compare this anisotropy
factor with the corresponding anisotropy factor from
ParFlow-input. Currently, we believe that the ParFlow-factor is used
in the beginning of the simulation, while the factor set here, is used
after each parameter update.

### PF:aniso_perm_z ###

`PF:aniso_perm_z`: (real) Anisotropy factor of saturated hydraulic
conductivity in z-direction. Only used when hydraulic conductivity is
updated (`[PF]paramupdate = 1`). Important: Compare this anisotropy
factor with the corresponding anisotropy factor from
ParFlow-input. Currently, we believe that the ParFlow-factor is used
in the beginning of the simulation, while the factor set here, is used
after each parameter update.

### PF:aniso_use_parflow ###

`PF:aniso_use_parflow`: (integer) If 1, heterogeneous anisotropy
factors (computed from ParFlow permeability arrays before the PDAF
update) are used to compute the permeability in y and z direction
based on the permeability in x direction that has been updated by
PDAF.

If 0, two homogeneous anisotropy factors are used to compute the
updated permeabilities in y and z direction. These factors are read
from the input variables `PF:aniso_perm_y` and `PF:aniso_perm_z`.

### PF:printensemble ###

`PF:printensemble`: (integer) If set to `1`, the updated state
variables for all ensemble members is printed out as `pfb` files after
each assimilation cycle.

Files are printed to `[DA]outdir` and follow the file naming
convention of ParFlow. They include the specifier `update` in the file
name.

For debug runs ([](precompiler:pdaf_debug)), the state ensemble is
also printed out before the DA update. The specifier for the files
containing the pre-DA state ensemble is `integrate`.

### PF:t_printensemble ###

`PF:t_printensemble`: (integer) The timestep for the state ensemble
output switched on under `PF:printensemble`.

Default setting is `-1`, which means: Print debug output at every DA
time step.

### PF:printstat ###

`PF:printstat`: (integer) If set to `1` (default) the ensemble
statistics (mean and standard deviation) of the forecasted state
variable are calculated and printed out as `pfb` files after each
assimilation cycle. Files are printed to `DA:outdir` and follow the
file naming convention of ParFlow. The output files include the
specifiers `press.mean`/ `press.sd` in case of assimilated pressure
update `PF:updateflag = 1` and the specifiers `swc.mean`/ `swc.sd` in
case of assimilated soil moisture data `PF:updateflag = 2` or
`PF:updateflag = 3`.

### PF:paramprintensemble ###

`PF:paramprintensemble`: (integer) Only used in case of parameter
update. If set to `1`, the updated parameters are printed to `pfb`
files (similar to `[PF]printensemble`). Output files include a
specifier depending on `PF:paramupdate`, f.e. `update.param.ksat`,
`update.param.ksat`.

### PF:paramprintstat ###

`PF:paramprintstat`: (integer) Only used in case of parameter
update. If set to `1` statistics on the updated parameters are printed
to `pfb` files (similar to `[PF]printstat`).

The output files include a specifier depending on `PF:paramupdate`,
f.e.  `param.ksat`, `param.mannings` or `param.poro`.

### PF:olfmasking ###

`PF:olfmasking`: (integer) Only used in case you do not want to update
the state on certain grid-cells during DA with pdaf. eg. not update
the cell which is saturated. Masked cells are not used for the state
update.

- Option \"1\": All cells at surface are masked.

- Option \"2\" reads a pfb of masked cells (use-case:
  rivers/streams). The full soil column below the masked cells is not
  updated!

- Option \"3\": All saturated cells at surface are masked. Only
  implemented, when pressure is in the state vector,
  i.e. `PF:updateflag` is `1` or `3`.

### PF:olfmasking_param ###

`PF:olfmasking_param`: (integer) Only used in case you do not want to
update the parameters from the EnKF state vector on certain grid-cells
during DA with pdaf. eg. not update the cell which is
saturated. Masked cells are not used for the parameter update.

- Option \"1\" All cells at surface are masked. Only implemented for
  `PF:paramupdate==1`.

- NOT YET IMPLEMENTED: Option \"2\" reads a pfb for masking the
  stream.

- Option \"3\": All saturated cells at surface are masked. Only
  implemented for `PF:paramupdate==1`; `PF:updateflag` `1` or `3`.

### PF:olfmasking_depth ###

`PF:olfmasking_depth`: (integer) Number of layers (counted from
surface) to mask in overland flow river masking. Default: 1.

(enkfpf:clm)=
## [CLM] ##

### CLM:problemname ###

`CLM:problemname`: (string) Problem prefix for CLM

(enkfpf:clm:nprocs)=
### CLM:nprocs ###

`CLM:nprocs`: (integer) Number of processors for each CLM instance.

This number of processors specifies a subset of the processors
available in a single `COMM_model`. Each `COMM_model` contains
`npes_world / nreal` processes.

### CLM:update_swc ###

`CLM:update_swc`: (integer) Flag for update of soil moisture content
in CLM (standalone only).

-  0: No update of soil moisture content

-  1: Update of soil moisture content

-  2: (only `clm3_5`) Update and average of soil moisture content for
   use with Cosmic-Ray data (Hui Pung implementation).

-  3: (only CLM5.0) Update and average of soil moisture content for
   use with Cosmic-Ray data (Strebel implementation of Schrön2017,
   <https://research-information.bris.ac.uk/en/publications/improving-calibration-and-validation-of-cosmic-ray-neutron-sensor>).

### CLM:update_texture ###

`CLM:update_texture`: (integer) Flag for update of soil parameter in
CLM (standalone only).

-  0: No update of soil parameter

-  1: Update of soil parameter (clay or sand)

-  (only CLM5.0) 2: Update of soil parameter (clay, sand or organic matter)

-  (only CLM5.0) 3: Update of hydraulic conductivity log-transformed

-  (only CLM5.0) 4: Update of hydraulic conductivity, porosity,
   suction, hydraulic conductivity exponent B (see CLM5.0 technical
   manual
   <https://escomp.github.io/ctsm-docs/versions/release-clm5.0/html/tech_note/index.html>)

### CLM:update_T ###

`CLM:update_T`: (integer) Flag for updating of ground and vegetation
temperature.

Currently only CLM3.5

-  0: No update of ground and vegetation temperature

-  1: Update of ground and vegetation temperature

### CLM:print_swc ###

`CLM:print_swc`: (integer) If set to `1`, the updated soil moisture
content in CLM for all ensemble members is printed out as `netcdf`
files (one file per realisation for the whole simulation
period). Files are printed to the run directory and are named
according to the CLM problem prefix of the corresponding realisation
and include the specifier `update` in the file name.

### CLM:print_et ###

`CLM:print_et`: (integer) Invoke function `write_clm_statistics`. For
further information, see source code. Default: `0`.

### CLM:statevec_allcol ###

`CLM:statevec_allcol`: (integer) Switch for using all SWC columns of a
eCLM gridcell in the state vector.

If `0` (default): Only one SWC value per grid cell is saved in the
state vector.

If `1`: `#columns` SWC values per grid cell are saved in the state
vector.

Remark: `CLM:statevec_allcol` should not be switched on together with
`CLM:statevec_colmean`.

### CLM:statevec_colmean ###

`CLM:statevec_colmean`: (integer) Switch for using gridcell averages
of SWC (which is a column-variable) in the state vector.

If `0` (default): Only the SWC value from the **first** column is
saved per grid cell in the state vector.

If `1`: The gridcell mean over all (possibly hydrologically active)
SWC values is saved per grid cell in the state vector.

Remark: `CLM:statevec_colmean` should not be switched on together with
`CLM:statevec_allcol`.

### CLM:statevec_only_active ###

`CLM:statevec_only_active`: (integer) Switch for using in the state
vector only (1) hydrologically active columns of a CLM5 gridcell and
(2) only layers until bedrock.

Only used for `CLM:update_swc`. Therein, only used, when either
`CLM:statevec_colmean` or `CLM:statevec_allcol` is used.

If `0` (default): Use all columns and all layers.

If `1`: Use only hydrologically active columns and only layers until
bedrock.

### CLM:statevec_max_layer ###

**Not yet in main branch**

`CLM:statevec_max_layer`: (integer) Number of layers to add in the
state vector.

Only used, when `CLM:statevec_allcol` and `CLM:statevec_only_active`
are switched on.

If `25` (default): All layers are in state vector.

If `9`: Only the first nine layers in state vector (corresponds to 1.2
meter).

For a depth profile of CLM layers, see [CLM Technical Note: 2.2.2.1
Soil
Layers](https://escomp.github.io/ctsm-docs/versions/master/html/tech_note/Ecosystem/CLM50_Tech_Note_Ecosystem.html#soil-layers).

### CLM:t_printensemble ###

`CLM:t_printensemble`: (integer) The timestep for the state ensemble
output switched on with the debug flag `PDAF_DEBUG`.

Default setting is `-1`, which means: Print debug output at every DA
time step.

### CLM:watmin_switch ###

`CLM:watmin_switch`: (integer) Switch for the values of minimal soil
moisture checked and, if lower, set during updating CLM's soil
moisture `h2osoi_vol`.

Default setting is `0`: Check and set `0.0`, i.e. no negative values
are allowed.

- `3`: CLM3.5 values: Check if SM in state vector is less than
  `0.00`. If yes, set SM to `0.05`.
- `5`: CLM5.0 values: Check if SM in state vector is less than
  CLM5.0's `watmin` from `clm_varcon.F90` (current value `0.01`). If
  yes, set SM to `watmin`.

(enkfpf:cosmo)=
## [COSMO] ##

(enkfpf:cosmo:nprocs)=
### COSMO:nprocs ###

`COSMO:nprocs`: (integer) Number of processors for each COSMO
instance.

Currently, `COSMO:nprocs` is NOT USED by TSMP-PDAF. Instead, COSMO
will get processors if there are processes left by ParFlow and CLM:
`(npes_world / nreal) - (PF:nprocs + CLM:nprocs)`.

### COSMO:dtmult ###

`COSMO:dtmult`: (integer) Number of COSMO time steps within one
ParFlow time step.

Remark: This holds for `PF:dt==1`, otherwise `PF:dt` is another factor
determining how COSMO rebuilds the ParFlow time step.

(enkfpf:da)=
## [DA] ##

### DA:outdir ###

`DA:outdir`: (string) Directory where assimilation results should be
written.

(enkfpf:da:nreal)=
### DA:nreal ###

`DA:nreal`: (integer) Number of realisations used in the
simulation. `DA:nreal` Must be equal to command line input
`n_modeltasks`.

Remark: Ensemble sizes of 2 have been shown to lead to
errors. Recommendation: Use at least an ensemble size of 4.

### DA:startreal ###

`DA:startreal`: (integer) Added to suffix-numbers for input file
creation.

### DA:total_steps ###

`DA:total_steps`: (integer) Number of observation Intervals.

Together with `DA:da_interval` and `DA:startcycle`, `DA:total_steps`
determines the simulation time.

If not set directly, `DA:total_steps` is determined from `PF:simtime`.

`DA:total_steps` must be in sync with timing information from
component models.

- ParFlow: `DA:total_steps` times `DA:da_interval` should be
  `TimingInfo.StopTime` minus `TimingInfo.StartTime` (from `*.pfidb`).

- CLM: `DA:total_steps` times `DA:da_interval` should be the number of
  CLM time steps in `stop_ymd` minus `start_ymd start_tod`. To be
  safe, CLM's stoptime can be chosen later, such that it will be
  stopped by the TSMP-PDAF stop alarm.

### DA:tstartcycle ###

`DA:tstartcycle`: (integer) First observation cycle (file) to use.

This should be zero. Other values are currently not supported.

(enkfpf:da:da_interval)=
### DA:da_interval ###

`DA:da_interval`: (double) Time interval (units of ParFlow timing,
usually hours, check ParFlow input `TimingInfo.BaseUnit`,
<https://parflow.readthedocs.io/en/latest/keys.html#timing-information>).

After the time interval `da_interval` the forward simulation
(integration) of all component models is stopped in the main loop of
TSMP-PDAF (`pdaf_terrsysmp.f90`). Afterwards PDAF is called and it
depends on the command line option `delt_obs` ([](cmd:delt_obs))
whether data assimilation is performed during this
iteration. `delt_obs` will specify the number of steps (loop
iterations) that PDAF will iterate back to the forward simulation,
before performing the actual data assimilation.

For CLM simulations, `DA:da_interval` determines the number of
CLM-time-steps between two assimilations together with `PF:dt`. See
the section of `PF:dt` for more detail. In particular, for CLM
simulations and `PF:dt==1.0`, `DA:da_interval` is in units of CLM time
steps.

One exception is that the observation file is empty at a data
assimilation step. Then, no data assimilation takes place and the next
set forward integration steps is executed.

Formula for the simulation time between data assimilation steps (given
non-empty observation files)

\begin{gather*} 
	t_{\mathtt{betweenDA}} =  \mathtt{da\_interval} \cdot \mathtt{delt\_obs} \cdot \mathtt{BaseUnit}
\end{gather*}


**Example 1 (FallSchoolCase):** `da_interval=1.0`, `delt_obs 12`,
`TimingInfo.BaseUnit=1.0`. In this case, component models will be
simulated for 12 1-hour-steps between data assimilation times. So an
assimilation will be applied each 12 hours.

**Example 2:** `da_interval=24.0`, `delt_obs 1`,
`TimingInfo.BaseUnit=1.0`. In this case, component models will be
simulated for 1 24-hour-step between data assimilation times. So an
assimilation will be applied each 24 hours.

**Remark**: Performance analysis has shown that for an assimilation
every `n` hours, it is beneficial to specify `da_interval=n`,
`delt_obs 1`, instead of `da_interval=1`, `delt_obs n`.

In general, it is beneficial to set `da_interval` as large as possible
for a given setup. One reason is that after each simulation time of
`da_interval`, the routines `assimilate_pdaf` and `update_tsmp` are
called, assembling EnKF state vectors and calling the PDAF
library. Maximizing `da_interval`, minimizes the number of these calls
and thus reduces compute time.

### DA:stat_dumpoffset ###

`DA:stat_dumpoffset`: File number offset for the data assimilation
output files for ParFlow.  Can be used when an assimilation run is
restarted.

### DA:da_interval_final ###

`DA:da_interval_final`: (float) Final `da_interval` for flexible time
stepping using observation files.

Used only if `DA:flexible_da_interval` is switched on.

Note: If the last observation file contains observations,
`da_interval_final` is set, but never used, as the simulation
currently finishes exactly at the time the observations from the last
observation file is assimilated.

### DA:flexible_da_interval ###

`DA:flexible_da_interval`: (integer) Switch for flexible time stepping
using observation files.

- `0`: default, fixed time stepping, single `da_interval`
- `1`: flexible time stepping, `da_interval` read from next
  observation file

- Set `total_steps` using `PF:simtime` (`DA:da_interval` must be `1`)
- Command line input [`delt_obs`](cmd:delt_obs) must be `1`!

In observation files, the new input variable
[`da_interval`](obs:timeflex:da_interval) determines the `da_interval`
before the observation file in question.

### DA:stat_dumpoffset ###

`DA:stat_dumpoffset`: File number offset for the data assimilation
output files for ParFlow.  Can be used when an assimilation run is
restarted.

(enkfpf:da:screen_wrapper)=
### DA:screen_wrapper ###

`DA:screen_wrapper` (added 02/2022): Variable for controlling
output. Analogous to `screen` variable for PDAF. Values: (0) no
outputs, (1) medium outputs (default), (2) debugging output.

(enkfpf:da:point_obs)=
### DA:point_obs ###

`DA:point_obs`: (integer) Switch to turn off point observations,
i.e. turn on [multi-scale data
assimilation](obs:files:clm:multi-scale-da).

Default value of `point_obs` is `1`. This means that point
observations are used for the data assimilation run.

If `point_obs` is set to `0`, multi-scale data assimilation is used.
Example: using SMAP satellite data over a large area.

(enkfpf:da:obs_interp_switch)=
### DA:obs_interp_switch ###

`DA:obs_interp_switch`: (integer) Switch for using an interpolation of
simulated measurements from the closest grid cells to the observation
location.

Default: `0` (no interpolation)

Effect of `obs_interp_switch=1`:
- For [CLM-type observations](obs:files:clm) the longitude and
  latitude inputs in the observation file are compared to the four
  neihboring grid points and the simulated measurements from these
  grid points are averaged with the distances to the grid points as
  weights.
- For [ParFlow-type observations](obs:files:pfl) the observation files
  need two more inputs: `ix_interp_d` and `iy_interp_d`. These index
  distances are used to compute an average of the simulated
  measurements from the value at the surrounding grid points in the
  x-y-plane.

### DA:crns_flag ###

`DA:crns_flag`: (int) Flag for using CRNS-observations

- `1`: ParFlow-based CRNS-measurement operator. ParFlow soil moisture
data is averaged and used in conventional weighting profile proposed
in Schroen etal HESS 2017 for modelling CRNS observations.

- `2`: CLM-based CRNS-measurement operator. CLM-soil moisture are
    weighted and used in the operator. CRNS implementation based on
    Schrön et al. 2017 using d86 for 3 different radius values. Note
    that the penetration depths are computed using the observations y.

Default: `0` (no CRNS-observations).

### DA:crns_bd ###

`DA:crns_bd`: (int) Only used for `DA:crns_flag=2`. Prescribed bulk
density average for the eight soil layers considered in
CRNS-measurement operator.

Default: `-1.0` (bulk density computed from CLM's
`soilstate_inst%bd_col`).

### DA:da_crns_depth_tol ###

`DA:da_crns_depth_tol`: (double) Only used for
`DA:crns_flag=1`. Convergence criteria for weighting procedure that
"generate[s] a weighted average of [soil-moisture] point measurements
that can be compared with [...] [a] CRNS product" (Schrön et al, 2017
http://dx.doi.org/10.5194/hess-21-5009-2017)

Only used if CRNS-observations are turned on by adding parameter
`depth` in observation input.

Relevant to assimilating PF averaged soil moisture data only for the
modelled CRNS data in Schroen etal in HESS 2017. 

This is a convergence criteria (between 0 and 1) for the averaging
skin depth (compare the weighting procedure described in Schrön et al,
2017, section 2.3).

### DA:print_obs_index ###

`DA:print_obs_index`: (int) Switch for turning on output of the
process-local grid-index-array of the obserations, `obs_index_p`.

For this output, debugging has to be turned on by `PDAF_DEBUG`.

Default: 0, output turned off.

(enkfpf:summary)=
## Parameter Summary ##

 | section   | parameter               | default value |
 |:---------:|:-----------------------:|:-------------:|
 | `[PF]`    |                         |               |
 |           | `problemname`           | \-            |
 |           | `nprocs`                | 0             |
 |           | `starttime`             | 0.0           |
 |           | `endtime`               | 0             |
 |           | `simtime`               | 0             |
 |           | `dt`                    | 0.0           |
 |           | `updateflag`            | 1             |
 |           | `paramupdate`           | 0             |
 |           | `paramupdate_frequency` | 1             |
 |           | `dampingfactor_param`   | 1.0           |
 |           | `dampingfactor_state`   | 1.0           |
 |           | `damping_switch_sm`     | 1             |
 |           | `aniso_perm_y`          | 1.0           |
 |           | `aniso_perm_z`          | 1.0           |
 |           | `aniso_use_parflow`     | 0             |
 |           |                         |               |
 |           | `printensemble`         | 1             |
 |           | `t_printensemble`       | -1            |
 |           | `printstat`             | 1             |
 |           | `paramprintensemble`    | 1             |
 |           | `paramprintstat`        | 1             |
 |           | `olfmasking`            | 0             |
 |           | `olfmasking_param`      | 0             |
 |           | `olfmasking_depth`      | 0             |
 | `[CLM]`   |                         |               |
 |           | `problemname`           | \-            |
 |           | `nprocs`                | 0             |
 |           | `update_swc`            | 1             |
 |           | `print_swc`             | 0             |
 |           | `print_et`              | 0             |
 |           | `statevec_allcol`       | 0             |
 |           | `statevec_colmean`      | 0             |
 |           | `statevec_only_active`  | 0             |
 |           | `statevec_max_layer`    | 25            |
 |           | `t_printensemble`       | -1            |
 |           | `watmin_switch`         | 0             |
 | `[COSMO]` |                         |               |
 |           | `nprocs`                | 0             |
 |           | `dtmult`                | 0             |
 | `[DA]`    |                         |               |
 |           | `nreal`                 | 0             |
 |           | `outdir`                | \-            |
 |           | `da_interval`           | 1             |
 |           | `da_interval_final`     | 1             |
 |           | `flexible_da_interval`  | 0             |
 |           | `stat_dumpoffset`       | 0             |
 |           | `screen_wrapper`        | 1             |
 |           | `point_obs`             | 1             |
 |           | `obs_interp_switch`     | 0             |
 |           | `crns_flag`             | 1             |
 |           | `da_crns_depth_tol`     | 0.01          |
 |           | `print_obs_index`       | 0             |

Default values for parameter file `enkfpf.par`.

