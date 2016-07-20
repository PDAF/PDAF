

/*-----------------------------------------------------------------------------------------
 * Copyright (c) 2013-2016 by Wolfgang Kurtz and Guowei He (Forschungszentrum Juelich GmbH)
 *
 * This file is part of TerrSysMP-PDAF
 *
 * TerrSysMP-PDAF is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * TerrSysMP-PDAF is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU LesserGeneral Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with TerrSysMP-PDAF.  If not, see <http://www.gnu.org/licenses/>.
 * -------------------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------------
 * enkf_parflow.h: Wrapper functions for ParFlow (header file)
 * -------------------------------------------------------------------------------------------*/


#ifndef GLOBAL
#define GLOBAL extern
#endif

#include "parflow.h"

#ifndef _ENKF_PARFLOW_H_
#define _ENKF_PARFLOW_H_

/* global integer variables */
GLOBAL int enkf_subvecsize;
GLOBAL int *idx_map_subvec2state;
GLOBAL int pf_statevecsize;
GLOBAL int pf_paramvecsize;
extern int pf_updateflag;
extern int pf_paramupdate;
GLOBAL int nx_local,ny_local,nz_local;
extern int pf_olfmasking;
GLOBAL int *riveridx,*riveridy,nriverid;

/* global double variables */
GLOBAL double *subvec_p, *subvec_sat, *subvec_porosity, * subvec_pressure_backup;
GLOBAL double *subvec_param;
GLOBAL double *pf_statevec;
GLOBAL double * xcoord, * ycoord, * zcoord;
extern double pf_aniso_perm_y,pf_aniso_perm_z;

/* global MPI communicator */
#if defined PARFLOW_STAND_ALONE
GLOBAL int comm_model_pdaf;
#endif

/* variables for calculation of statistics */
GLOBAL double *subvec_mean, *subvec_sd;
extern int    comm_couple, task_id;

/* functions */
void enkfparflowinit(int ac, char *av[],char *input_file); 
void enkfparflowadvance(double current_time, double dt);
void enkfparflowfinalize();
void enkf_printvec(char *pre, char *suff, double *data);
void enkf_printmannings(char *pre, char *suff);
void enkf_ensemblestatistics (double* dat, double* mean, double* var, int size, MPI_Comm comm);
void parflow_oasis_init(double current_time, double dt);
void init_idx_map_subvec2state(Vector *pf_vector);

void PF2ENKF(Vector *pf_vector, double *enkf_subvec);
void ENKF2PF(Vector *pf_vector, double *enkf_subvec);
int  enkf_getsubvectorsize(Grid *grid);

void update_parflow();
void mask_overlandcells();
void mask_overlandcells_river();

#endif


