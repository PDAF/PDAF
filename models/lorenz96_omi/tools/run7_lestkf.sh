#!/bin/sh
# This script is an example showing how to run a sequence of assimilation with
# varying parameters.
# In this example we run the SEIK filter with ensemble size 30. We perform a
# set of experiments with varying values for the forgetting factor.
#
# 2010-02, Lars Nerger, AWI
# $Id: runasml.sh 350 2020-01-30 18:35:17Z lnerger $

# Name of executable
EXE="./pdaf_lorenz_96"

# General settings for all experiments
DEFAULTS="-total_steps 2 -step_null 1000 -dim_ens 8"

FILTER=7

# Run experiments
for FORGET in 0.95 #1 0.99 0.98 0.97 0.96 0.95 0.94 0.93 0.92 0.91 0.9
do
    $EXE $DEFAULTS -filtertype $FILTER -forget $FORGET -local_range 10 -locweight 2 \
	-file_asml t${FILTER}_N30_f${FORGET}.nc
done

