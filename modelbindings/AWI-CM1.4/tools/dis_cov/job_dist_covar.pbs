#!/bin/bash
#PBS -N dist_covar
#PBS -j oe
#PBS -l walltime=01:00:00
#PBS -l nodes=1:ppn=1
#PBS -q mpp2testq
##PBS -l feature=ice2 
#PBS -A hbk00064

#module load mvapich2/1.2.0-intel
#module load netcdf
#module load intel.compiler
#module load cray-parallel-netcdf/1.8.1.3

set -vx
export F_UFMTENDIAN=big
ulimit -s

cd $PBS_O_WORKDIR

date
aprun -n 1 ./distribute_covar_AWI-CM
date

