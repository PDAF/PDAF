#!/bin/bash
#PBS -N eofcovar
#PBS -j oe
#PBS -l walltime=00:10:00
#PBS -l nodes=1:ppn=1
#PBS -q mpp2testq 
##PBS -l feature=ice2 
##PBS -A hbkqtang
#PBS -A hbk00064

#module load mvapich2/1.2.0-intel
#module load netcdf
module load cray-parallel-netcdf/1.8.1.3
#module load intel.compiler
module load intel/18.0.0.128

set -vx
export F_UFMTENDIAN=big

cd $PBS_O_WORKDIR

date
aprun -n 1 ./generate_covar_AWI-CM
date

