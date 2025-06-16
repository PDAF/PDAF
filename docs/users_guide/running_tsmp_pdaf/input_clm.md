(clm)=
# Input files for CLM #

For executing TSMP-PDAF, the standard CLM input file `lnd.stdin` is
replaced by separate input files for the individual CLM realisations.
The naming convention for these input files is a user defined prefix
followed by a formatted realisation number. For example, if the chosen
prefix is `clminput`, the file names for the first three realisations
are:

```text
clminput_00000
clminput_00001
clminput_00002
```

All these files follow the same structure as the standard CLM input
file `lnd.stdin`. The timing information needs to be the same for all
these CLM input files but all other variables representing CLM input
and output can differ between the realisations. It is important to
note that the variable `caseid` should be an individual identifier for
each CLM realisation because this variable determines the names of the
CLM output files. If certain realisation share the same `caseid` the
name of the CLM output files will be the same for these realisation
(i.e., these realisation will produce the same output files) which
results in undefined behaviour.

