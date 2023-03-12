#!/bin/sh
# This script is an example to run the foreward Lorenz96 model
# for 10000 time steps.
#
# 2010-02, Lars Nerger, AWI
# $Id: runmodel.sh 79 2019-02-09 15:29:27Z lnerger $

# Name of executable
EXE="./lorenz_63"

# General settings for all experiments
DEFAULTS="-total_steps 10000"

# Run model
$EXE $DEFAULTS
