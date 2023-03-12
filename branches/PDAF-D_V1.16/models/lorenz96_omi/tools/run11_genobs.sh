#!/bin/sh
# This script is an example showing how we can generate
# synthetic observations with PDAF.
#
# We run with ensemble size 2 here to avoid that NaNs
# appear in the output. However, running with dim_ens=1
# also works and is recommended.
#
# 2019-02, Lars Nerger, AWI
# $Id: rungenobs.sh 180 2019-07-03 15:19:46Z lnerger $

# Name of executable
EXE="./pdaf_lorenz_96"

# General settings for all experiments
DEFAULTS="-total_steps 20 -step_null 1000 -dim_ens 2 -type_ensinit 'tru'"

# Run experiments
FORGET=1 
$EXE $DEFAULTS -filtertype 11 -forget $FORGET \
	-file_asml t11_N2_genobs.nc
