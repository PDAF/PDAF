#!/bin/bash
#$Id: runtests_offline_crayCS400.slurm 1703 2016-12-17 12:02:33Z lnerger $
#############################################################
### Script for testing different configuration options    ###
### of PDAF. This version is to test all subtypes of the  ###
### different filters with parallelization.               ###
### The particular configuration is for the Cray CS400    ###
### system at AWI with slurm batch system.                ###
###                                                       ###
### The data assimilation is run with 2 model tasks,      ###
### where each model uses domain decomposition. The cases ###
### with a static covariance matrix use only 1 task.      ###
#############################################################

# Settings for batch system
#SBATCH --job-name=pdaftest
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p mpp
#SBATCH --time=00:10:00

export OMP_NUM_THREADS=1

# General configuration
NENS=50                    # Ensemble size in EnKF/SEIK/LSEIK
NEOF=`expr $NENS - 1`      # Number of EOFs in SEEK
CONF="-dim_state 300 -screen 1 -tasks 1"          # General configuration for dynamic filters
CONF_FIXED="-dim_state 300 -screen 1 -tasks 1"    # General configuration for fixed covariance
EXE="./pdaf_dummy_offline"  # Name of executable
CMD="srun"                 # Command for parallel execution

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
$CMD $EXE $CONF -dim_ens $NEOF -filtertype 0 -subtype 5 -filename output_par_seek5.dat
fi
if [ $TEST_SEIK -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run SEIK tests"
echo "--------------------------------------------------------"
$CMD $EXE $CONF -dim_ens $NENS -filtertype 1 -subtype 5 -filename output_par_seik5.dat
fi
if [ $TEST_ENKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run EnKF tests"
echo "--------------------------------------------------------"
$CMD $EXE $CONF -dim_ens $NENS -filtertype 2 -subtype 5 -filename output_par_enkf5.dat
fi
if [ $TEST_LSEIK -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run LSEIK tests"
echo "--------------------------------------------------------"
$CMD $EXE $CONF -dim_ens $NENS -filtertype 3 -subtype 5 -filename output_par_lseik5.dat
fi
if [ $TEST_ETKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run ETKF tests"
echo "--------------------------------------------------------"
$CMD $EXE $CONF -dim_ens $NENS -filtertype 4 -subtype 5 -filename output_par_etkf5.dat
fi
if [ $TEST_LETKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run LETKF tests"
echo "--------------------------------------------------------"
$CMD $EXE $CONF -dim_ens $NENS -filtertype 5 -subtype 5 -filename output_par_letkf5.dat
fi
if [ $TEST_ESTKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run ESTKF tests"
echo "--------------------------------------------------------"
$CMD $EXE $CONF -dim_ens $NENS -filtertype 6 -subtype 5 -filename output_par_estkf5.dat
fi
if [ $TEST_LESTKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run LESTKF tests"
echo "--------------------------------------------------------"
$CMD $EXE $CONF -dim_ens $NENS -filtertype 7 -subtype 5 -filename output_par_lestkf5.dat
fi
if [ $TEST_LENKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run LENKF tests"
echo "--------------------------------------------------------"
$CMD $EXE $CONF -dim_ens $NENS -filtertype 8 -subtype 5 -filename output_par_lenkf5.dat
fi
if [ $TEST_NETF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run NETF tests"
echo "--------------------------------------------------------"
$CMD $EXE $CONF -dim_ens $NENS -filtertype 9 -subtype 5 -filename output_par_netf5.dat
fi
if [ $TEST_LNETF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run LNETF tests"
echo "--------------------------------------------------------"
$CMD $EXE $CONF -dim_ens $NENS -filtertype 10 -subtype 5 -filename output_par_lnetf5.dat
fi

echo "PDAF tests completed: " `date`


exit
