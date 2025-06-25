(cmd)=
# Command line options #

The following command line options must be specified when TSMP-PDAF is
executed:

## n_modeltasks ##

`n_modeltasks` (integer) Number of realisations. Must be equal to
[`[DA] nreal`](enkfpf:da:nreal).

Note that both `n_modeltasks` (or `nreal`) and the total number of MPI
tasks specified, f.e. with `srun -n` or `mpiexec -np` (internally:
`npes_world`) have to match with the total number of processes
specified for the component models:
[`PF:nprocs`](enkfpf:pf:nprocs),
[`CLM:nprocs`](enkfpf:clm:nprocs),
[`COSMO:nprocs`](enkfpf:cosmo:nprocs).

## filtertype ##

`filtertype` (integer) Type of filter used for data assimilation. For
more details see the PDAF documentation. Currently, the following
filter options are available

-   (Ensemble Kalman Filter(Enkf)) is supported.

-   (Ensemble Transform Kalman Filter(Etkf)) is supported.

-   (Local Ensemble Transform Kalman Filter(Letkf)) is supported.

-   (Error Subspace Transform Kalman Filter(Estkf)) is supported.

-   (Local Error Subspace Transform Kalman Filter(Lestkf)) is supported.

-   (Local Ensemble Kalman Filter(LEnkf)) is supported.

To see all available command line options for the different filter
algorithms please refer to PDAF wiki link

<http://pdaf.awi.de/trac/wiki/AvailableOptionsforInitPDAF>

## subtype ##

`subtype` (integer) Parameter subtype, different options for each
filter. See [](cmd:command-line-examples).

## obs_filename ##

`obs_filename` (string) Prefix for observation files.

Note: The prefix for observation files cannot be longer than 100
characters! Otherwise, the observation filename will be cropped.

Tip: Try to use relative paths as input to `obs_filename`.

If observation filenames are not found, an error message `No such file
or directory` is printed out in `next_observation_pdaf.F90`.
 
## rms_obs ##

`rms_obs` (real) Measurement error. This is the standard deviation of
the Gaussian probability distribution of the observations.

`rms_obs` should be given in units of the observation.

(cmd:delt_obs)=
## delt_obs ##

`delt_obs` (integer) Number of data assimilation intervals (see
[`[DA]da_interval`](enkfpf:da:da_interval)) that are forward computed,
before assimilation is actually performed.

In general, `delt_obs` should be as small as possible, in order to
avoid performance loss. See remark in
[`[DA]da_interval`](enkfpf:da:da_interval).

For flexible time stepping, `delt_obs` must be one.

(cmd:screen)=
## screen ##

`screen` (integer) Control verbosity of PDAF
- 0: no outputs
- 1: basic output (default)
- 2: 1 plus timing output
- 3: 2 plus debug output

## forget ##

`forget` (real) forgetting factor for filter analysis

Example: `-forget 0.98`.

General advise: Choose forgetting factor close to one. For values
smaller than 0.95, effects like a splitting of the ensemble have been
observed (compare Amezcua et al., Tellus A 2012, 64, 18039,
<http://dx.doi.org/10.3402/tellusa.v64i0.18039>)

For EnKF / LEnKF, the forgetting factor leads to a spreading of the
ensemble through manipulating ensemble member by

\begin{align*}
x^{f}_{i} &= \bar{x} + (x_{i}-\bar{x}) \cdot \frac{1}{\mathtt{forget}^2},
\end{align*}

where $x_{i}$ is the state vector ensemble member $i$ and $\bar{x}$ is
the ensemble mean of the state vector.

For ETKF, see
e.g. <https://github.com/PDAF/PDAF/blob/ae9545227bd4804469dff389a9baadcc9e31906e/src/PDAF_etkf_analysis.F90#L441-L444>

## locweight ##

Only for localization.

Usage depends on interface routines.

`locweight` (integer): set weight function for localization,
- default=0
- for constant weight of 1
- possible are integer values 0 to 4.(see `init_pdaf`)

``` text
locweight == 0    ! Uniform (unit) weighting
locweight == 1    ! Exponential weighting
locweight == 2    ! 5th-order polynomial (Gaspari&Cohn, 1999)
```
## cradius ##

Only for localization.

Usage depends on interface routines.

`cradius` (old name: `local_range`) (real) set cut-off radius
- 0.0 by default
- any positive value should work.

## sradius ##

Only for localization.

Usage depends on interface routines.

`sradius` (old name `srange`) (real): the support radius of the
localization, default: equal to `cradius`.

`sradius` has different meanings depending on the value of
`locweight`. 

- `locweight == 0`: `sradius` is not used
- `locweight == 1`: `sradius` is distance with weight `1/e`, for
  avoiding sharp cut-off, `cradius` should be significantly larger
  than `sradius`
- `locweight == 2`: `sradius` is the support radius for 5th-order
  polynomial, can safely be chosen the same as `cradius`,
  parametrization change of the 5th-order polynomial is at `sradius/2`,
  see
  <https://github.com/PDAF/PDAF/blob/ed631034956dece8e91e8b588c4cf3aaa7916f49/src/PDAF_local_weight.F90#L147-L176>.

More detail in PDAF documentation:
<https://pdaf.awi.de/trac/wiki/PDAF_local_weight>


(cmd:command-line-examples)=
## Command Line Examples ##

### Command Line Example: EnKF ###

Main option: `-filtertype 2`

Additional options:
- `-subtype`

Example of the model execution using EnKF analysis:

        mpiexec -np 512 ./tsmp-pdaf -n_modeltasks 64 -filtertype 2 -subytype 1 -delt_obs 1 -rms_obs 0.1 -obs_filename myobs

With this command, the TSMP-PDAF executable `tsmp-pdaf` will be run with
64 realisations (8 processors for each realisation) and EnKF analysis
will be performed every assimilation cycle with an observation error of
0.1 and observations read from the files `myobs.xxxxx`.

#### EnKF: subtype ####

`subtype` (integer) Parameter subtype for Ensemble Kalman Filter (EnKF).

``` fortran
     IF (subtype == 0) THEN
        WRITE (*, '(a, 14x, a)') 'PDAF', '--> EnKF with analysis for large observation dimension'
     ELSE IF (subtype == 1) THEN
        WRITE (*, '(a, 14x, a)') 'PDAF', '--> EnKF with analysis for small observation dimension'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(a, 14x, a)') 'PDAF', '--> offline mode'

```


### Command Line Example: ETKF ###

Main option: `-filtertype 4`

Additional options:
- `-subtype`

Example of the model execution using EtKF analysis:

        mpiexec -np 512 ./tsmp-pdaf -n_modeltasks 64 -filtertype 4 -subytype 1 -delt_obs 1 -rms_obs 0.1 -obs_filename myobs

With this command, the TSMP-PDAF executable `tsmp-pdaf` will be run with
64 realisations (8 processors for each realisation) and EtKF analysis
will be performed every assimilation cycle with an observation error of
0.1 and observations read from the files `myobs.xxxxx`.

#### ETKF: subtype ####

`subtype` (integer) Parameter subtype for Ensemble Transform Kalman
Filter (ETKF).

``` fortran
     IF (subtype == 0) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> ETKF using T-matrix'
     ELSE IF (subtype == 1) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> ETKF following Hunt et al. (2007)'
     ELSE IF (subtype == 2) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> ETKF with fixed error-space basis'
     ELSE IF (subtype == 3) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> ETKF with fixed state covariance matrix'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> offline mode'

```

### Command Line Example: ESTKF ###

Main option: `-filtertype 6`

Additional options:
- `-subtype`

Example of the model execution using EstKF analysis:

        mpiexec -np 512 ./tsmp-pdaf -n_modeltasks 64 -filtertype 6 -subytype 1 -delt_obs 1 -rms_obs 0.1 -obs_filename myobs

With this command, the TSMP-PDAF executable `tsmp-pdaf` will be run with
64 realisations (8 processors for each realisation) and EstKF analysis
will be performed every assimilation cycle with an observation error of
0.1 and observations read from the files `myobs.xxxxx`.

#### ESTKF: subtype ####

`subtype` (integer) Parameter subtype for Error Subspace Transform Kalman Filter
(ESTKF).

``` fortran
     IF (subtype == 0) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> Standard ESTKF'
     ELSE IF (subtype == 2) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> ESTKF with fixed error-space basis'
     ELSE IF (subtype == 3) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> ESTKF with fixed state covariance matrix'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> offline mode'

```

### Command Line Example: LETKF ###

Main option: `-filtertype 5`

Additional options:
- `-subtype`
- `-locweight`
- `-cradius`
- `-sradius`

Example of the model execution using Letkf analysis:

        mpiexec -np 512 ./tsmp-pdaf -n_modeltasks 64 -filtertype 5 -cradius 3 -locweight 0 -subtype 0 -delt_obs 1 -rms_obs 0.1 -obs_filename myobs

With this command, the TSMP-PDAF executable `tsmp-pdaf` will be run with
64 realisations (8 processors for each realisation) and Letkf analysis
will be performed every assimilation cycle with an observation error of
0.1 and observations read from the files `myobs.xxxxx`.

#### LETKF subtype

`subtype` (integer) Parameter subtype for Local Ensemble Transform
Kalman Filter (LETKF).

```f90
     IF (subtype == 0) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> LETKF using T-matrix'
     ELSE IF (subtype == 1) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> LETKF following Hunt et al. (2007)'
     ELSE IF (subtype == 2) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> LETKF with fixed error-space basis'
     ELSE IF (subtype == 3) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> LETKF with fixed state covariance matrix'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> offline mode'
```

### Command Line Example: LESTKF ###

Main option: `-filtertype 7`

Additional options:
- `-subtype`
- `-locweight`
- `-cradius`
- `-sradius`

Example of the model execution using LESTKF analysis:

        mpiexec -np 512 ./tsmp-pdaf -n_modeltasks 64 -filtertype 7 -cradius 3 -locweight 0 -subtype 0 -delt_obs 1 -rms_obs 0.1 -obs_filename myobs

With this command, the TSMP-PDAF executable `tsmp-pdaf` will be run with
64 realisations (8 processors for each realisation) and Lestkf analysis
will be performed every assimilation cycle with an observation error of
0.1 and observations read from the files `myobs.xxxxx`.

#### LESTKF: subtype ####

`subtype` (integer) Parameter subtype for Local Error Subspace Transform Kalman
Filter (LESTKF).


``` fortran
     IF (subtype == 0) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> Standard LESTKF'
     ELSE IF (subtype == 2) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> LESTKF with fixed error-space basis'
     ELSE IF (subtype == 3) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> LESTKF with fixed state covariance matrix'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> offline mode'
```

### Command Line Example: LENKF ###

Main option: `-filtertype 8`

Additional options:
- `-subtype`
- `-locweight`
- `-cradius`
- `-sradius`

**Example** of the model execution using Lenkf analysis:

``` bash
mpiexec -np 512 ./tsmp-pdaf -n_modeltasks 64 -filtertype 8 -subtype 0 -locweight 2 -cradius 3 -delt_obs 1 -rms_obs 0.1 -obs_filename myobs
```

With this command, the TSMP-PDAF executable `tsmp-pdaf` will be run
with 64 realisations (8 processors for each realisation) and Lenkf
analysis will be performed every assimilation cycle with an
observation error of 0.1 and observations read from the files
`myobs.xxxxx`.

#### LENKF: subtype ####

`subtype` (integer) Parameter subtype for Localized Ensemble Kalman
Filter (Lenkf).
-   Full ensemble integration; analysis with covariance localization
-   Offline mode
