#!/bin/sh
# In this example we run the ESTKF filter with ensemble size 20 and
# a forgetting factor of 0.8. The run is a twin experiment using the
# synthetic observation which were generated before using rungenobs.sh
#
# 2019-02, Lars Nerger, AWI
# $Id: runtwin.sh 79 2019-02-09 15:29:27Z lnerger $

# Name of executable
EXE="./pdaf_lorenz_63"

# General settings for all experiments
DEFAULTS="-total_steps 5000 -dim_ens 20 -forget 0.8"

# Run experiments
$EXE $DEFAULTS -filtertype 6 -file_asml ESTKF_N20_twin.nc -twin_experiment t


