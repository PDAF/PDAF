(cos)=
# Input files for COSMO #

Currently, different COSMO realizations share the same input namelist
files (`INPUT_`), except the namelist file for input/output
specifications (`INPUT_IO`). The `INPUT_IO` namelist file names for
the first three model realisation are:

``` text
INPUT_IO_00000
INPUT_IO_00001
INPUT_IO_00002
```

Within these different IO namelist files it is possible to assign
different directories for the model input and output of the respective
model realization.

