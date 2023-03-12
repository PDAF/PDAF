#!/bin/sh
#$Id: runtests_offline_1pe.sh 1700 2016-12-17 11:01:23Z lnerger $

#############################################################
### Script for testing different configuration options    ###
### of PDAF. This version is to test the offline mode for ###
### the different filters without parallelization.        ###
#############################################################

#set -vx

# General configuration
NENS=50                 # Ensemble size in EnKF/SEIK/LSEIK
NEOF=`expr $NENS - 1`    # Number of EOFs in SEEK
CONF="-dim_state 300 -screen 1"    # General configuration (state dimension)
EXE="./pdaf_dummy_offline"     # Name of executable

TEST_SEEK=1   # (1) to perform tests with the SEEK filter
TEST_SEIK=1   # (1) to perform tests with the SEIK filter
TEST_ENKF=1   # (1) to perform tests with the Ensemble Kalman filter
TEST_LSEIK=1  # (1) to perform tests with the local SEIK filte
TEST_ETKF=1   # (1) to perform tests with the ETKF
TEST_LETKF=1  # (1) to perform tests with the LETKF
TEST_ESTKF=1  # (1) to perform tests with the ESTKF
TEST_LESTKF=1 # (1) to perform tests with the LESTKF
TEST_LENKF=1  # (1) to perform tests with the localized Ensemble Kalman filter
TEST_NETF=1   # (1) to perform tests with the NETF
TEST_LNETF=1  # (1) to perform tests with the LNETF

# Perform tests
echo "================  Testing PDAF offline ================="

echo "Machine: " `uname -a`
echo "Date: " `date`

if [ $TEST_SEEK -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run SEEK tests"
echo "--------------------------------------------------------"
$EXE $CONF -dim_ens $NEOF -filtertype 0 -subtype 5 -filename output_seek5.dat
fi
if [ $TEST_SEIK -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run SEIK tests"
echo "--------------------------------------------------------"
$EXE $CONF -dim_ens $NENS -filtertype 1 -subtype 5 -filename output_seik5.dat
fi
if [ $TEST_ENKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run EnKF tests"
echo "--------------------------------------------------------"
$EXE $CONF -dim_ens $NENS -filtertype 2 -subtype 5 -filename output_enkf5.dat
fi
if [ $TEST_LSEIK -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run LSEIK tests"
echo "--------------------------------------------------------"
$EXE $CONF -dim_ens $NENS -filtertype 3 -subtype 5 -filename output_lseik5.dat
fi
if [ $TEST_ETKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run ETKF tests"
echo "--------------------------------------------------------"
$EXE $CONF -dim_ens $NENS -filtertype 4 -subtype 5 -filename output_etkf5.dat
fi
if [ $TEST_LETKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run LETKF tests"
echo "--------------------------------------------------------"
$EXE $CONF -dim_ens $NENS -filtertype 5 -subtype 5 -filename output_letkf5.dat
fi
if [ $TEST_ESTKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run ESTKF tests"
echo "--------------------------------------------------------"
$EXE $CONF -dim_ens $NENS -filtertype 6 -subtype 5 -filename output_estkf5.dat
fi
if [ $TEST_LESTKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run LESTKF tests"
echo "--------------------------------------------------------"
$EXE $CONF -dim_ens $NENS -filtertype 6 -subtype 5 -filename output_lestkf5.dat
fi
if [ $TEST_LENKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run LEnKF tests"
echo "--------------------------------------------------------"
$EXE $CONF -dim_ens $NENS -filtertype 8 -subtype 5 -filename output_lenkf5.dat
fi
if [ $TEST_NETF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run NETF tests"
echo "--------------------------------------------------------"
$EXE $CONF -dim_ens $NENS -filtertype 9 -subtype 5 -filename output_netf5.dat
fi
if [ $TEST_LENKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run LNETF tests"
echo "--------------------------------------------------------"
$EXE $CONF -dim_ens $NENS -filtertype 10 -subtype 5 -filename output_lnetf5.dat
fi

echo "PDAF tests completed: " `date`

