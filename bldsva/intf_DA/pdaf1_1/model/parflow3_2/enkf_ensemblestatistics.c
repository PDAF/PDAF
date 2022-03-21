/*-----------------------------------------------------------------------------------------
Copyright (c) 2013-2016 by Wolfgang Kurtz and Guowei He (Forschungszentrum Juelich GmbH)

This file is part of TerrSysMP-PDAF

TerrSysMP-PDAF is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

TerrSysMP-PDAF is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU LesserGeneral Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with TerrSysMP-PDAF.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------------
enkf_ensemblestatistics.c: Functions for calculating ensemble statistics for ParFlow
-------------------------------------------------------------------------------------------*/

#include "enkf.h"
#include "enkf_parflow.h"
#include <math.h>

/*-------------------------------------------------------------------------*/
/**
  @author   Wolfgang Kurtz, Guowei He
  @brief    Compute ensemble statistics
  @param[in]    double* dat Input vector (ensembles distributed across `comm`)
  @param[out]   double* mean Mean vector of ensemble `dat` computed in this routine
  @param[out]   double* var Variance vector of ensemble `dat` computed in this routine
  @param[in]    int size Size of `mean`, `var` and a single realization `dat`
  @param[in]    MPI_Comm comm Communicator for MPI-reduction-operation

  1. Sum up ensemble members from different PEs (`MPI_Allreduce`)
  2. Normalize sums to obtain mean; compute variance summands
  3. Sum up variance summands (`MPI_Reduce`)
  4. Normalize variance sums to obtain variance
 */
/*--------------------------------------------------------------------------*/
void enkf_ensemblestatistics (double* dat, double* mean, double* var, int size, MPI_Comm comm)
{
  int i;
  double *varsum;

  varsum = (double*) calloc(size,sizeof(double));

  /* 1. Sum of ensemble members */
  MPI_Allreduce(dat,mean,size,MPI_DOUBLE,MPI_SUM,comm);

  for(i=0;i<size;i++){
      /* 2. Normalize sum to obtain mean */
      mean[i] = mean[i] / nreal;
      /* 2. Compute variance summands */
      var[i]  = (dat[i] - mean[i]) * (dat[i] - mean[i]);
  }

  //MPI_Allreduce(MPI_IN_PLACE,var,size,MPI_DOUBLE,MPI_SUM,comm);
  /* 3. Sum up variance summands */
  MPI_Reduce(var,varsum,size,MPI_DOUBLE,MPI_SUM,0,comm);

  /* 4. Normalize variance sums to obtain variance */
  for(i=0;i<size;i++) var[i] = sqrt(varsum[i] / (nreal-1));

  free(varsum);
}

void enkf_printstatistics_pfb (double *dat, char* name, int cycle, char* prefix, int dim)
{
  char outfile[200];
  char outfile_ts[10];

  //sprintf(outfile,"%s/%s.%s",dir,pfinfile,name);
  sprintf(outfile,"%s.%s",prefix,name);
  sprintf(outfile_ts,"%05d",cycle);
  enkf_printvec(outfile,outfile_ts,dat,dim);
}
