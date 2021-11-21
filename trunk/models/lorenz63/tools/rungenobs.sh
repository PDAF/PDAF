#!/bin/sh
# This script is an example showing how we can generate
# synthetic observations with PDAF.
#
# We run with ensemble size 2 here to avoid that NaNs
# appear in the output. However, running with dim_ens=1
# also works and is recommended.
#
# 2019-02, Lars Nerger, AWI
# $Id: rungenobs.sh 79 2019-02-09 15:29:27Z lnerger $

# Name of executable
EXE="./pdaf_lorenz_63"

# General settings for all experiments
DEFAULTS="-total_steps 10000 -step_null 1000 -dim_ens 2 -type_ensinit 'tru'"

# Run experiments
FORGET=1 
$EXE $DEFAULTS -filtertype 11 -file_asml t11_N1_genobs.nc -use_obs_mask f
