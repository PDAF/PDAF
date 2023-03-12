#!/bin/sh
# This script is an example to run the foreward Lorenz96 model
# for 10000 time steps.
#
# 2010-02, Lars Nerger, AWI
# $Id: runasml.sh 860 2010-02-03 10:05:15Z lnerger $

# Name of executable
EXE="./lorenz_96"

# General settings for all experiments
DEFAULTS="-total_steps 10000"

# Run model
$EXE $DEFAULTS
