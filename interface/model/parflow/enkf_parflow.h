/*-----------------------------------------------------------------------------------------
Copyright (c) 2013-2016 by Wolfgang Kurtz, Guowei He and Mukund Pondkule (Forschungszentrum Juelich GmbH)

This file is part of TSMP-PDAF

TSMP-PDAF is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

TSMP-PDAF is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU LesserGeneral Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with TSMP-PDAF.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------------
enkf_parflow.h: Wrapper functions for ParFlow (header file)
-------------------------------------------------------------------------------------------*/


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
GLOBAL int nx_glob,ny_glob,nz_glob;
GLOBAL int nx_local,ny_local,nz_local;
GLOBAL int nx_glob, ny_glob, nz_glob;
GLOBAL int origin_local[3];
extern int pf_olfmasking;
extern int pf_olfmasking_param;
extern int pf_olfmasking_depth;
extern int pf_gwmasking;
extern int pf_printgwmask;
GLOBAL int *riveridx,*riveridy,nriverid;

/* global double variables */
GLOBAL double *subvec_p, *subvec_sat, *subvec_porosity, *subvec_param;
GLOBAL double *subvec_gwind;
GLOBAL double *pf_statevec;
GLOBAL double *subvec_permy;
GLOBAL double *subvec_permz;
GLOBAL double *arr_aniso_perm_yy;
GLOBAL double *arr_aniso_perm_zz;
GLOBAL double * xcoord;
GLOBAL double * ycoord;
GLOBAL double * zcoord;
/* hcp CRNS begins */
//GLOBAL double *subvec_Kind;    //hcp
GLOBAL double *soilay;  //hcp soil layers
/* hcp CRNS ends */
extern double pf_aniso_perm_y,pf_aniso_perm_z;

/* global MPI communicator */
#if defined PARFLOW_STAND_ALONE
GLOBAL int COMM_model_pfl;
#endif

/* variables for calculation of statistics */
GLOBAL double *subvec_mean, *subvec_sd;
GLOBAL double *subvec_param_mean, *subvec_param_sd;
extern int    comm_couple;  /* task_id; */
GLOBAL double *dat_alpha, *dat_n, *dat_ksat, *dat_poro;

/* functions */
void enkfparflowinit(int ac, char *av[],char *input_file);
void enkfparflowadvance(int tcycle, double current_time, double dt);
void enkfparflowfinalize();
void enkf_printvec(char *pre, char *suff, double *data, int dim);
void enkf_printmannings(char *pre, char *suff);
void enkf_ensemblestatistics (double* dat, double* mean, double* var, int size, MPI_Comm comm);
void parflow_oasis_init(double current_time, double dt);
void init_idx_map_subvec2state(Vector *pf_vector);

void PF2ENKF(Vector *pf_vector, double *enkf_subvec);
void PF2ENKF_2P(Vector *p1_vector, Vector *p2_vector, double *enkf_subvec);
void PF2ENKF_3P(Vector *p1_vector, Vector *p2_vector, Vector *p3_vector, double *enkf_subvec);
void PF2ENKF_4P(Vector *p1_vector, Vector *p2_vector, Vector *p3_vector, Vector *p4_vector, double *enkf_subvec);
void ENKF2PF(Vector *pf_vector, double *enkf_subvec);
void ENKF2PF_2P(Vector *p1_vector, Vector *p2_vector, double *enkf_subvec);
void ENKF2PF_3P(Vector *p1_vector, Vector *p2_vector, Vector *p3_vector, double *enkf_subvec);
void ENKF2PF_4P(Vector *p1_vector, Vector *p2_vector, Vector *p3_vector, Vector *p4_vector, double *enkf_subvec);
void ENKF2PF_masked(Vector *pf_vector, double *enkf_subvec, double *mask);
int  enkf_getsubvectorsize(Grid *grid);

void update_parflow();
void mask_overlandcells();
void mask_overlandcells_river();
void init_n_domains_size(int* n_domains_p);
void init_parf_l_size(int* dim_l);
//void g2l_state(int* domain_p, float* state_p[], int* dim_l, float* state_l[]);
//void l2g_state(int* domain_p, float* state_p[], int* dim_l, float* state_l[]);
void print_update_pfb();

/* external functions/ variables (fortran/ pdaf) for retrieving measurement locations for current time step */
extern void get_obsindex_currentobsfile(int *no_obs);
extern void clean_obs_pf();
extern int *tidx_obs, *xidx_obs, *yidx_obs, *zidx_obs, *ind_obs;

#endif
