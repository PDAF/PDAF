#!/bin/sh
# In this example we run the ESTKF filter with ensemble size 20 and
# a forgetting factor of 0.8.
#
# 2019-07, Lars Nerger, AWI
# $Id: runasml.sh 79 2019-02-09 15:29:27Z lnerger $

# Name of executable
EXE="./pdaf_lorenz_63"

# General settings for all experiments
DEFAULTS="-total_steps 5000 -dim_ens 20 -forget 0.8"

# Run experiments
$EXE $DEFAULTS -filtertype 6 -file_asml ESTKF_N20.nc


