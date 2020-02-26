#!/bin/sh
# This script is an example showing how to run a sequence of assimilation with
# varying parameters.
# In this example we run the particle filter with ensemble size 20.
# Further we use residual resampling (pf_noise_type 2) and
# ensemble perturbations with noise of amplitude 0.2 (pf_noise_amp 0.2)
# relative to the ensmeble variance (pf_noise_type 2)
#
# 2019-07, Lars Nerger, AWI
# $Id: runasml.sh 79 2019-02-09 15:29:27Z lnerger $

# Name of executable
EXE="./pdaf_lorenz_63"

# General settings for all experiments
DEFAULTS="-total_steps 5000 -dim_ens 20 -pf_res_type 2 -pf_noise_type 2 -pf_noise_amp 0.2"

# Run experiments
$EXE $DEFAULTS -filtertype 12 -file_asml PF_N20.nc


