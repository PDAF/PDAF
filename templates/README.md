Collection of template routines for PDAF user routines

## Directories

offline/ 
- offline coupling using PDAF3 interface

online/
- online coupling using PDAF3 interface in fully parallel mode

online_flexible/
- online coupling using PDAF3 interface in flexible parallel mode

omi/
- templates specific for PDAF-OMI (file with callback routines, observation module, observation operator)

3dvar/
- additional templates for implementing 3D-Var methods (for online implementation; first copy online/, then copy 3dvar/ overwriting some files)

full_interface/
- classical PDAF-1 implementation using the interface routines specifying all call-back routines (not using OMI or PDAFlocal)


## Compiling and running

One can compile the online and offline tutorial cases in each directory like the tutorial cases. Note that while compiling and running is possible, there is no useful functionality in the user-supplied call-back routines. However, one can base ones own implementation on the templates and do a step-wise implementation as described in the tutorials. For this, the templates include screen outputs that remind you which routines are not yet completed.

The online cases can be run e.g. with
    mpirun -np 4 ./PDAF_online -dim_ens 4 
which runs the ESTKF. In addition one could specify the filter type with '-filtertype X' where e.g. X=6 or X=7 can be specified.

The offline template can be run with
    ./PDAF_offline -dim_ens 4 
or simply by
	./PDAF_offline 
in which case the ESTKF with dim_ens=9 is used from the template code. In addition, one could specify the filter type with '-filtertype X' where e.g. X=6 or X=7 can be specified.
