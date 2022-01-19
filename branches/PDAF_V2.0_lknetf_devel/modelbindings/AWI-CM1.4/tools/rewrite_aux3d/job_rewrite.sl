#!/bin/bash 
##SBATCH --account=hbk00064
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
##SBATCH --output=20ens-%j.out
##SBATCH --error=20ens-%j.err
#SBATCH --time=00:30:00
#SBATCH --partition=standard96:test
##SBATCH --mail-user=lars.nerger@awi.de
##SBATCH --mail-type=END
#SBATCH -A zzz0002

module load intel/18.0.6
module load impi

#set -vx
export F_UFMTENDIAN=big
ulimit -s

export SLURM_CPU_BIND=none
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/sw/dataformats/netcdf/intel.18/4.7.3/skl/lib

date
mpirun -np 1 ./rewrite_aux3d
date

