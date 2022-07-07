#!/bin/sh
# This script is an example showing how to run a sequence of assimilation with
# varying parameters.
# In this example we run the SEIK filter with ensemble size 30. We perform a
# single experiment a forgetting facrtor of 0.95.
#
# 2021-11, Lars Nerger, AWI
# $Id: runtwin.sh 350 2020-01-30 18:35:17Z lnerger $

# Name of executable
EXE="./pdaf_lorenz_05c"

# General settings for all experiments
DEFAULTS="-total_steps 500 -step_null 1000 -dim_ens 30"

# Run experiments
for FORGET in 0.95 
do
    $EXE $DEFAULTS -filtertype 6 -forget $FORGET -twin_experiment t \
	-file_asml t1_N30_f${FORGET}_twin.nc
done

