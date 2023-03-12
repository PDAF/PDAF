#!/usr/bin/sh
#$Id:$
#############################################################
### Script for testing different configuration options    ###
### of PDAF. This version is to test the offline mode for ###
### the different filters without parallelization.        ###
### The particular configuration is for the IBM Regatta   ###
### at AWI with load leveler batch system.                ###
#############################################################

# Setings for batch system
# @ job_name = pdaftest
# @ step_name = $(job_name)
# @ output = out.runtests_offline_ibm.$(jobid)
# @ error = $(output)
# @ step_name = $(job_name)
# @ notification = never
# @ job_type = parallel
# @ class = cpar
# @ node_usage = not_shared
# @ node = 1
# @ tasks_per_node = 4
# @ resources = ConsumableCpus(1)
# @ wall_clock_limit = 3600
# @ queue

# General configuration
NENS=50                  # Ensemble size in EnKF/SEIK/LSEIK
NEOF=`expr $NENS - 1`    # Number of EOFs in SEEK
CONF="-dim_state 300 -screen 1"    # General configuration (state dimension)
EXE="./pdaf_dummy_offline"         # Name of executable
CMD="poe"                # Command for parallel execution

TEST_SEEK=1   # (1) to perform tests with the SEEK filter
TEST_SEIK=1   # (1) to perform tests with the SEIK filter
TEST_ENKF=1   # (1) to perform tests with the Ensemble Kalman filter
TEST_LSEIK=1  # (1) to perform tests with the local SEIK filter
TEST_ETKF=1   # (1) to perform tests with the ETKF
TEST_LETKF=1  # (1) to perform tests with the LETKF
TEST_ESTKF=1  # (1) to perform tests with the ESTKF
TEST_LESTKF=1 # (1) to perform tests with the LESTKF


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

echo "PDAF tests completed: " `date`


exit
