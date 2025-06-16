# Input files for OASIS-MCT #

When the component models of TSMP are coupled via OASIS-MCT, certain
input files for the coupler need to be provided for a regular TSMP run.
These files are principally the same for the execution of TSMP-PDAF and
are shared among the different realisations. These input files include
the file `namcouple` which defines the coupling between CLM and ParFlow,
the netCDF file `clmgrid.nc` which is identical to the CLM input file
containing the grid data (see CLM documentation) and the remapping files
for the coupling (`rmp_.nc`). See the OASIS-MCT documentation for more
details. It is important to note here that the remapping files should be
created by a single deterministic TSMP run of the model. The so created
remapping files should then be copied to the TSMP-PDAF run directory. If
the remapping files are not present when TSMP-PDAF is executed, each
realisation will try to create these files. However, because the naming
convention for the OASIS-MCT remapping files is fixed, this would lead
to an undefined overwriting of these files from the different
realisations.

