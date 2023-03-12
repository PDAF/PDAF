#!/bin/sh
#$Id: runtests_1pe.sh 1700 2016-12-17 11:01:23Z lnerger $

#############################################################
### Script for testing different configuration options    ###
### of PDAF. This version is to test all subtypes of the  ###
### different filters without parallelization.            ###
#############################################################

#set -vx

# General configuration
NENS=50                 # Ensemble size in EnKF/SEIK/LSEIK
NEOF=`expr $NENS - 1`    # Number of EOFs in SEEK
CONF="-dim_state 70 -screen 1"    # General configuration (state dimension)
EXE="./pdaf_dummy_snglprec"  # Name of executable

TEST_SEEK=1   # (1) to perform tests with the SEEK filter
TEST_SEIK=1   # (1) to perform tests with the SEIK filter
TEST_ENKF=1   # (1) to perform tests with the Ensemble Kalman filter
TEST_LSEIK=1  # (1) to perform tests with the local SEIK filter
TEST_ETKF=1   # (1) to perform tests with the ETKF
TEST_LETKF=1  # (1) to perform tests with the LETKF
TEST_ESTKF=1  # (1) to perform tests with the ESTKF
TEST_LESTKF=1 # (1) to perform tests with the LESTKF
TEST_LENKF=1  # (1) to perform tests with the localized Ensemble Kalman filter
TEST_NETF=1   # (1) to perform tests with the NETF
TEST_LNETF=1  # (1) to perform tests with the LNETF

# Perform tests
echo "====================  Testing PDAF  ===================="

echo "Machine: " `uname -a`
echo "Date: " `date`

if [ $TEST_SEEK -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run SEEK tests"
echo "--------------------------------------------------------"
$EXE $CONF -dim_ens $NEOF -filtertype 0 -subtype 0 -filename output_seek0.dat
$EXE $CONF -dim_ens $NEOF -filtertype 0 -subtype 1 -filename output_seek1.dat
$EXE $CONF -dim_ens $NEOF -filtertype 0 -subtype 2 -filename output_seek2.dat
$EXE $CONF -dim_ens $NEOF -filtertype 0 -subtype 3 -filename output_seek3.dat
fi
if [ $TEST_SEIK -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run SEIK tests"
echo "--------------------------------------------------------"
$EXE $CONF -dim_ens $NENS -filtertype 1 -subtype 0 -filename output_seik0.dat
$EXE $CONF -dim_ens $NENS -filtertype 1 -subtype 1 -filename output_seik1.dat
$EXE $CONF -dim_ens $NENS -filtertype 1 -subtype 2 -filename output_seik2.dat
$EXE $CONF -dim_ens $NENS -filtertype 1 -subtype 3 -filename output_seik3.dat
$EXE $CONF -dim_ens $NENS -filtertype 1 -subtype 4 -filename output_seik4.dat
fi
if [ $TEST_ENKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run EnKF tests"
echo "--------------------------------------------------------"
$EXE $CONF -dim_ens $NENS -filtertype 2 -subtype 0 -filename output_enkf0.dat
$EXE $CONF -dim_ens $NENS -filtertype 2 -subtype 1 -filename output_enkf1.dat
fi
if [ $TEST_LSEIK -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run LSEIK tests"
echo "--------------------------------------------------------"
$EXE $CONF -dim_ens $NENS -filtertype 3 -subtype 0 -filename output_lseik0.dat
$EXE $CONF -dim_ens $NENS -filtertype 3 -subtype 2 -filename output_lseik2.dat
$EXE $CONF -dim_ens $NENS -filtertype 3 -subtype 3 -filename output_lseik3.dat
$EXE $CONF -dim_ens $NENS -filtertype 3 -subtype 4 -filename output_lseik4.dat
fi
if [ $TEST_ETKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run ETKF tests"
echo "--------------------------------------------------------"
$EXE $CONF -dim_ens $NENS -filtertype 4 -subtype 0 -filename output_etkf0.dat
$EXE $CONF -dim_ens $NENS -filtertype 4 -subtype 1 -filename output_etkf1.dat
$EXE $CONF -dim_ens $NENS -filtertype 4 -subtype 2 -filename output_etkf2.dat
$EXE $CONF -dim_ens $NENS -filtertype 4 -subtype 3 -filename output_etkf3.dat
fi
if [ $TEST_LETKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run LETKF tests"
echo "--------------------------------------------------------"
$EXE $CONF -dim_ens $NENS -filtertype 5 -subtype 0 -filename output_letkf0.dat
$EXE $CONF -dim_ens $NENS -filtertype 5 -subtype 1 -filename output_letkf1.dat
$EXE $CONF -dim_ens $NENS -filtertype 5 -subtype 2 -filename output_letkf2.dat
$EXE $CONF -dim_ens $NENS -filtertype 5 -subtype 3 -filename output_letkf3.dat
fi
if [ $TEST_ESTKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run ESTKF tests"
echo "--------------------------------------------------------"
$EXE $CONF -dim_ens $NENS -filtertype 6 -subtype 0 -filename output_estkf0.dat
$EXE $CONF -dim_ens $NENS -filtertype 6 -subtype 2 -filename output_estkf2.dat
$EXE $CONF -dim_ens $NENS -filtertype 6 -subtype 3 -filename output_estkf3.dat
fi
if [ $TEST_LESTKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run LESTKF tests"
echo "--------------------------------------------------------"
$EXE $CONF -dim_ens $NENS -filtertype 7 -subtype 0 -filename output_lestkf0.dat
$EXE $CONF -dim_ens $NENS -filtertype 7 -subtype 2 -filename output_lestkf2.dat
$EXE $CONF -dim_ens $NENS -filtertype 7 -subtype 3 -filename output_lestkf3.dat
fi
if [ $TEST_LENKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run LEnKF tests"
echo "--------------------------------------------------------"
$EXE $CONF -dim_ens $NENS -filtertype 8 -subtype 0 -filename output_lenkf0.dat
fi
if [ $TEST_NETF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run NETF tests"
echo "--------------------------------------------------------"
$EXE $CONF -dim_ens $NENS -filtertype 9 -subtype 0 -filename output_netf0.dat
fi
if [ $TEST_LNETF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run LNETF tests"
echo "--------------------------------------------------------"
$EXE $CONF -dim_ens $NENS -filtertype 10 -subtype 0 -filename output_lnetf0.dat
fi

echo "PDAF tests completed: " `date`

