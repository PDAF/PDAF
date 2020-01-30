#!/bin/sh
# This script is an example showing how to run a sequence of assimilation with
# varying parameters.
# In this example we run the SEIK filter with ensemble size 30. We perform a
# single experiment a forgetting facrtor of 0.95.
#
# 2019-02, Lars Nerger, AWI
# $Id$

# Name of executable
EXE="./pdaf_lorenz_96"

# General settings for all experiments
DEFAULTS="-total_steps 500 -step_null 1000 -dim_ens 30"

# Run experiments
for FORGET in 0.95 
do
    $EXE $DEFAULTS -filtertype 6 -forget $FORGET -twin_experiment t \
	-file_asml t1_N30_f${FORGET}_twin.nc
done

