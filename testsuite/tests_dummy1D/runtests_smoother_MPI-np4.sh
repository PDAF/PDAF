#!/bin/sh
#$Id: runtests_smoother_MPI-np4.sh 1700 2016-12-17 11:01:23Z lnerger $

#############################################################
### Script for testing different configuration options    ###
### of PDAF. This version is to test all subtypes of the  ###
### different smoothers with parallelization.             ###
###                                                       ###
### The data assimilation is run with 2 model tasks,      ###
### where each model uses domain decomposition with 2     ###
### processes per model.                                  ###
#############################################################

#set -vx

# General configuration
NENS=50                    # Ensemble size in EnKF/SEIK/LSEIK
NEOF=`expr $NENS - 1`      # Number of EOFs in SEEK
CONF="-dim_state 300 -screen 1 -tasks 2 -dim_lag 10"    # General configuration (state dimension)
CONF_FIXED="-dim_state 300 -screen 1 -tasks 1 -dim_lag 10"    # General configuration for fixed covariance
EXE="./pdaf_dummy_online"  # Name of executable
CMD="mpirun -np 4"         # Command for parallel execution

TEST_ENKF=1   # (1) to perform tests with the Ensemble Kalman smooother
TEST_ETKF=1   # (1) to perform tests with the ETKF smoother
TEST_LETKF=1  # (1) to perform tests with the LETKF smoother
TEST_ESTKF=1  # (1) to perform tests with the ESTKF smoother
TEST_LESTKF=1 # (1) to perform tests with the LESTKF smoother
TEST_NETF=1   # (1) to perform tests with the NETF smoother
TEST_LNETF=1  # (1) to perform tests with the LNETF smoother

# Perform tests
echo "====================  Testing smoother in PDAF  ===================="

echo "Machine: " `uname -a`
echo "Date: " `date`

if [ $TEST_ENKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run EnKF tests"
echo "--------------------------------------------------------"
$CMD $EXE $CONF -dim_ens $NENS -filtertype 2 -subtype 0 -filename output_par_smoother_enkf0.dat
fi
if [ $TEST_ETKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run ETKF tests"
echo "--------------------------------------------------------"
$CMD $EXE $CONF -dim_ens $NENS -filtertype 4 -subtype 0 -filename output_par_smoother_etkf0.dat
$CMD $EXE $CONF -dim_ens $NENS -filtertype 4 -subtype 1 -filename output_par_smoother_etkf1.dat
$CMD $EXE $CONF_FIXED -dim_ens $NENS -filtertype 4 -subtype 2 -filename output_par_smoother_etkf2.dat
fi
if [ $TEST_LETKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run LETKF tests"
echo "--------------------------------------------------------"
$CMD $EXE $CONF -dim_ens $NENS -filtertype 5 -subtype 0 -filename output_par_smoother_letkf0.dat
$CMD $EXE $CONF -dim_ens $NENS -filtertype 5 -subtype 1 -filename output_par_smoother_letkf1.dat
$CMD $EXE $CONF_FIXED -dim_ens $NENS -filtertype 5 -subtype 2 -filename output_par_smoother_letkf2.dat
fi
if [ $TEST_ESTKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run ESTKF tests"
echo "--------------------------------------------------------"
$CMD $EXE $CONF -dim_ens $NENS -filtertype 6 -subtype 0 -filename output_par_smoother_estkf0.dat
$CMD $EXE $CONF_FIXED -dim_ens $NENS -filtertype 6 -subtype 2 -filename output_par_smoother_estkf2.dat
fi
if [ $TEST_LESTKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run LESTKF tests"
echo "--------------------------------------------------------"
$CMD $EXE $CONF -dim_ens $NENS -filtertype 7 -subtype 0 -filename output_par_smoother_lestkf0.dat
$CMD $EXE $CONF_FIXED -dim_ens $NENS -filtertype 7 -subtype 2 -filename output_par_smoother_lestkf2.dat
fi
if [ $TEST_NETF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run NETF tests"
echo "--------------------------------------------------------"
$CMD $EXE $CONF -dim_ens $NENS -filtertype 9 -subtype 0 -filename output_par_smoother_netf0.dat
fi
if [ $TEST_LNETF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run LNETF tests"
echo "--------------------------------------------------------"
$CMD $EXE $CONF -dim_ens $NENS -filtertype 10 -subtype 0 -filename output_par_smoother_netf0.dat
fi

echo "PDAF tests completed: " `date`

