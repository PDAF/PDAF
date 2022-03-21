/*-----------------------------------------------------------------------------------------
Copyright (c) 2013-2016 by Wolfgang Kurtz, Guowei He and Mukund Pondkule (Forschungszentrum Juelich GmbH)

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
enkf.h: Header file for global variables/ functions
-------------------------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>

#ifndef GLOBAL
#define GLOBAL extern
#endif

/* functions */
void read_enkfpar(char *parname);
void enkf_ensemblestatistics (double* dat, double* mean, double* var, int size, MPI_Comm comm);
void enkf_printstatistics_pfb (double *dat, char* name, int cycle, char* prefix, int dim);
extern void clm_init(char *s);
extern void clm_advance(int *ntstep);
extern void clm_finalize();
extern void cosmo_init();
extern void cosmo_advance(int *cos_dt);
extern void cosmo_finalize();

/* chars */
GLOBAL char pfinfile[500];
GLOBAL char pfprefixin[100];
GLOBAL char pfprefixout[100];
GLOBAL char pfoutfile_ens[500];
GLOBAL char pfoutfile_stat[500];
GLOBAL char pfproblemname[100];
GLOBAL char clminfile[100*2];
GLOBAL char outdir[100];

/* integers */
GLOBAL int nprocpf;
GLOBAL int nprocclm;
GLOBAL int nproccosmo;
GLOBAL int nreal;
GLOBAL int startreal;
GLOBAL int total_steps;
GLOBAL int stat_dumpint;
GLOBAL int stat_dumpoffset;
GLOBAL int screen_wrapper;
GLOBAL int point_obs;
GLOBAL MPI_Fint fsubcomm;
GLOBAL int oasprefixno;
GLOBAL int clmprefixlen;
GLOBAL int pf_updateflag;
GLOBAL int pf_paramupdate;
GLOBAL int pf_printensemble;
GLOBAL int pf_printstat;
GLOBAL int pf_paramprintensemble;
GLOBAL int pf_paramprintstat;
GLOBAL int nx_local,ny_local,nz_local;
GLOBAL int clmupdate_swc;
GLOBAL int clmupdate_T;
GLOBAL int clmupdate_texture;
GLOBAL int clmprint_swc;
GLOBAL int clmprint_T;   // hcp
GLOBAL int clmprint_et;
GLOBAL int dtmult_cosmo;
GLOBAL int pf_olfmasking;
GLOBAL int pf_gwmasking;
GLOBAL int pf_printgwmask;
GLOBAL int pf_freq_paramupdate;
extern int model;
extern int mype_model;
extern int npes_model;
extern int mype_world;
extern int npes_world;
extern int mype_filter;
extern int npes_filter;
extern int task_id;
extern int n_modeltasks;

/* double */
GLOBAL double *pmean,*satmean,*pvar,*satvar;
GLOBAL double t_start,t_end,dt;
GLOBAL double pf_aniso_perm_y,pf_aniso_perm_z;
GLOBAL double da_interval;
GLOBAL double pf_dampfac_param;
