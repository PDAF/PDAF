#!/bin/sh
# This script is an example showing how to run a sequence of assimilation with
# varying parameters.
# In this example we run the SEIK filter with ensemble size 30. We perform a
# set of experiments with varying values for the forgetting factor.
#
# 2010-02, Lars Nerger, AWI
# $Id: runasml.sh 350 2020-01-30 18:35:17Z lnerger $

./tools/run1_seik.sh
./tools/run2_enkf.sh
./tools/run3_lseik.sh
./tools/run4_etkf.sh
./tools/run5_letkf.sh
./tools/run6_estkf.sh
./tools/run7_lestkf.sh
./tools/run8_lenkf.sh
./tools/run9_netf.sh
./tools/run10_lnetf.sh
./tools/run11_genobs.sh
./tools/run12_pf.sh
