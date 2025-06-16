# Running TSMP-PDAF #

The execution of TSMP-PDAF and the necessary input is closely related
to the one of TSMP. For executing TSMP, a ParFlow database file
`*.pfidb` (see ParFlow documentation for more details), an input file
for CLM (`lnd.stdin`) and several input files for OASIS-MCT and COSMO
need to be present in the run directory. TSMP-PDAF follows a similar
approach. The only difference is that for each realisation a different
set of ParFlow/ CLM/ COSMO input files needs to be present which
follow a certain naming convention (see Sections [Input files for
CLM](clm), [Input files for ParFlow](pfl) and [Input files for
COSMO](cos)). Additionally, a control file for the data assimilation
([Control file `enkfpf.par`](enkfpf)) and observation files
([Observation files](obs)) need to be present in the run directory.
Furthermore, some command line options ([Command line options](cmd))
need to be specified when TSMP-PDAF is executed.

See the Virtual Machine download on webpage
<https://datapub.fz-juelich.de/slts/tsmp-vm/index.html>.

The data on this website is quite large - `12.8 GB`.

Only download this, if you really want to run the Virtual Machine!.


