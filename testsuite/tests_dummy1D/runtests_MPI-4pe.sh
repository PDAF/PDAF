#!/bin/bash
#$Id: runtests_MPI-4pe.sh 1700 2016-12-17 11:01:23Z lnerger $
#############################################################
### Script for testing different configuration options    ###
### of PDAF. This version is to test all subtypes of the  ###
### different filters with parallelization.               ###
### This is a general configuration to test with MPI.     ###
###                                                       ###
### The data assimilation is run with 2 model tasks,      ###
### where each model uses domain decomposition with 2     ###
### processes per model. The cases with a static          ###
### covariance matrix use only 1 task.                    ###
#############################################################

# General configuration
NENS=50                    # Ensemble size in EnKF/SEIK/LSEIK
NEOF=`expr $NENS - 1`      # Number of EOFs in SEEK
CONF="-dim_state 300 -screen 1 -tasks 2"          # General configuration for dynamic filters
CONF_FIXED="-dim_state 300 -screen 1 -tasks 1"    # General configuration for fixed covariance
EXE="./pdaf_dummy_online"  # Name of executable
CMD="mpirun -np 4"         # Command for parallel execution

TEST_SEEK=1   # (1) to perform tests with the SEEK filter
TEST_SEIK=1   # (1) to perform tests with the SEIK filter
TEST_ENKF=1   # (1) to perform tests with the Ensemble Kalman filter
TEST_LSEIK=1  # (1) to perform tests with the local SEIK filter
TEST_ETKF=1   # (1) to perform tests with the ETKF
TEST_LETKF=1  # (1) to perform tests with the LETKF
TEST_ESTKF=1  # (1) to perform tests with the ESTKF
TEST_LESTKF=1 # (1) to perform tests with the LESTKF
TEST_LENKF=1  # (1) to perform tests with the localized EnKF
TEST_NETF=1   # (1) to perform tests with the NETF
TEST_LNETF=1  # (1) to perform tests with the localized NETF

# Perform tests
echo "====================  Testing PDAF  ===================="

echo "Machine: " `uname -a`
echo "Date: " `date`

if [ $TEST_SEEK -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run SEEK tests"
echo "--------------------------------------------------------"
$CMD $EXE $CONF -dim_ens $NEOF -filtertype 0 -subtype 0 -filename output_par_seek0.dat
$CMD $EXE $CONF -dim_ens $NEOF -filtertype 0 -subtype 1 -filename output_par_seek1.dat
$CMD $EXE $CONF_FIXED -dim_ens $NEOF -filtertype 0 -subtype 2 -filename output_par_seek2.dat
$CMD $EXE $CONF_FIXED -dim_ens $NEOF -filtertype 0 -subtype 3 -filename output_par_seek3.dat
fi
if [ $TEST_SEIK -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run SEIK tests"
echo "--------------------------------------------------------"
$CMD $EXE $CONF -dim_ens $NENS -filtertype 1 -subtype 0 -filename output_par_seik0.dat
$CMD $EXE $CONF -dim_ens $NENS -filtertype 1 -subtype 1 -filename output_par_seik1.dat
$CMD $EXE $CONF_FIXED -dim_ens $NENS -filtertype 1 -subtype 2 -filename output_par_seik2.dat
$CMD $EXE $CONF_FIXED -dim_ens $NENS -filtertype 1 -subtype 3 -filename output_par_seik3.dat
$CMD $EXE $CONF -dim_ens $NENS -filtertype 1 -subtype 4 -filename output_par_seik4.dat
fi
if [ $TEST_ENKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run EnKF tests"
echo "--------------------------------------------------------"
$CMD $EXE $CONF -dim_ens $NENS -filtertype 2 -subtype 0 -filename output_par_enkf0.dat
$CMD $EXE $CONF -dim_ens $NENS -filtertype 2 -subtype 1 -filename output_par_enkf1.dat
fi
if [ $TEST_LSEIK -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run LSEIK tests"
echo "--------------------------------------------------------"
$CMD $EXE $CONF -dim_ens $NENS -filtertype 3 -subtype 0 -filename output_par_lseik0.dat
$CMD $EXE $CONF_FIXED -dim_ens $NENS -filtertype 3 -subtype 2 -filename output_par_lseik2.dat
$CMD $EXE $CONF_FIXED -dim_ens $NENS -filtertype 3 -subtype 3 -filename output_par_lseik3.dat
$CMD $EXE $CONF -dim_ens $NENS -filtertype 3 -subtype 4 -filename output_par_lseik4.dat
fi
if [ $TEST_ETKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run ETKF tests"
echo "--------------------------------------------------------"
$CMD $EXE $CONF -dim_ens $NENS -filtertype 4 -subtype 0 -filename output_par_etkf0.dat
$CMD $EXE $CONF -dim_ens $NENS -filtertype 4 -subtype 1 -filename output_par_etkf1.dat
$CMD $EXE $CONF_FIXED -dim_ens $NENS -filtertype 4 -subtype 2 -filename output_par_etkf2.dat
$CMD $EXE $CONF_FIXED -dim_ens $NENS -filtertype 4 -subtype 3 -filename output_par_etkf3.dat
fi
if [ $TEST_LETKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run LETKF tests"
echo "--------------------------------------------------------"
$CMD $EXE $CONF -dim_ens $NENS -filtertype 5 -subtype 0 -filename output_par_letkf0.dat
$CMD $EXE $CONF -dim_ens $NENS -filtertype 5 -subtype 1 -filename output_par_letkf1.dat
$CMD $EXE $CONF_FIXED -dim_ens $NENS -filtertype 5 -subtype 2 -filename output_par_letkf2.dat
$CMD $EXE $CONF_FIXED -dim_ens $NENS -filtertype 5 -subtype 3 -filename output_par_letkf3.dat
fi
if [ $TEST_ESTKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run ESTKF tests"
echo "--------------------------------------------------------"
$CMD $EXE $CONF -dim_ens $NENS -filtertype 6 -subtype 0 -filename output_par_estkf0.dat
$CMD $EXE $CONF_FIXED -dim_ens $NENS -filtertype 6 -subtype 2 -filename output_par_estkf2.dat
$CMD $EXE $CONF_FIXED -dim_ens $NENS -filtertype 6 -subtype 3 -filename output_par_estkf3.dat
fi
if [ $TEST_LESTKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run LESTKF tests"
echo "--------------------------------------------------------"
$CMD $EXE $CONF -dim_ens $NENS -filtertype 7 -subtype 0 -filename output_par_lestkf0.dat
$CMD $EXE $CONF_FIXED -dim_ens $NENS -filtertype 7 -subtype 2 -filename output_par_lestkf2.dat
$CMD $EXE $CONF_FIXED -dim_ens $NENS -filtertype 7 -subtype 3 -filename output_par_lestkf3.dat
fi
if [ $TEST_LENKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run LENKF tests"
echo "--------------------------------------------------------"
$CMD $EXE $CONF -dim_ens $NENS -filtertype 8 -subtype 0 -filename output_par_lenkf0.dat
fi
if [ $TEST_NETF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run NETF tests"
echo "--------------------------------------------------------"
$CMD $EXE $CONF -dim_ens $NENS -filtertype 9 -subtype 0 -filename output_par_netf0.dat
fi
if [ $TEST_LNETF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run NETF tests"
echo "--------------------------------------------------------"
$CMD $EXE $CONF -dim_ens $NENS -filtertype 10 -subtype 0 -filename output_par_lnetf0.dat
fi

echo "PDAF tests completed: " `date`


exit
