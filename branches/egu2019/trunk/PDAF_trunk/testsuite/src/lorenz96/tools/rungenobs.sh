#!/bin/sh
# This script is an exam[le showing how one can generate synthetic observations
# with PDAF. 
#
# 2019-02, Lars Nerger, AWI
# $Id$

# Name of executable
EXE="./pdaf_lorenz_96"

# General settings for all experiments
DEFAULTS="-total_steps 2000 -step_null 1000 -dim_ens 1 -type_ensinit 'tru'"

# Run experiments
FORGET=1 
$EXE $DEFAULTS -filtertype 11 -forget $FORGET \
	-file_asml t11_N1_genobs.nc
