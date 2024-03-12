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
read_enkfpar.c: Function for reading controle file of TSMP-PDAF
-------------------------------------------------------------------------------------------*/

#include "enkf.h"
#include "iniparser.h"

void read_enkfpar(char *parname)
{
  char *string;
  dictionary *pardict;

  int coupcol;
 
  /* initialize dictionary */
  pardict = iniparser_load(parname);
 
  /* get settings for ParFlow */
  string                = iniparser_getstring(pardict,"PF:problemname", "");
  strcpy(pfproblemname,string);
  string                = iniparser_getstring(pardict,"PF:prefix_input","");
  strcpy(pfprefixin,string);
  string                = iniparser_getstring(pardict,"PF:prefix_output","");
  strcpy(pfprefixout,string);
  nprocpf               = iniparser_getint(pardict,"PF:nprocs",0);
  t_start               = iniparser_getdouble(pardict,"PF:starttime",0);
  t_sim                 = iniparser_getdouble(pardict,"PF:simtime",0);
  dt                    = iniparser_getdouble(pardict,"PF:dt",0);
  pf_updateflag         = iniparser_getint(pardict,"PF:updateflag",1);
  pf_paramupdate        = iniparser_getint(pardict,"PF:paramupdate",0);
  pf_aniso_perm_y       = iniparser_getdouble(pardict,"PF:aniso_perm_y",1);
  pf_aniso_perm_z       = iniparser_getdouble(pardict,"PF:aniso_perm_z",1);
  pf_aniso_use_parflow  = iniparser_getint(pardict,"PF:aniso_use_parflow",0);
  pf_printensemble      = iniparser_getint(pardict,"PF:printensemble",1);
  pf_printstat          = iniparser_getint(pardict,"PF:printstat",1);
  pf_paramprintensemble = iniparser_getint(pardict,"PF:paramprintensemble",1);
  pf_paramprintstat     = iniparser_getint(pardict,"PF:paramprintstat",1);
  pf_olfmasking         = iniparser_getint(pardict,"PF:olfmasking",0);
  pf_olfmasking_param   = iniparser_getint(pardict,"PF:olfmasking_param",0);
  pf_olfmasking_depth   = iniparser_getint(pardict,"PF:olfmasking_depth",1);
  pf_gwmasking          = iniparser_getint(pardict,"PF:gwmasking",0);
  pf_printgwmask        = iniparser_getint(pardict,"PF:printgwmask",0);
  pf_dampfac_param      = iniparser_getdouble(pardict,"PF:dampingfactor_param",1.0);
  pf_dampfac_state      = iniparser_getdouble(pardict,"PF:dampingfactor_state",1.0);
  pf_dampswitch_sm        = iniparser_getdouble(pardict,"PF:damping_switch_sm",0);
  pf_freq_paramupdate   = iniparser_getint(pardict,"PF:paramupdate_frequency",1);

  /* backward compatibility settings for ParFlow */
  if (t_sim == 0){
    t_sim                 = iniparser_getdouble(pardict,"PF:endtime",0);
  }
  
  /* get settings for CLM */
  string                = iniparser_getstring(pardict,"CLM:problemname", "");
  strcpy(clminfile,string);
  nprocclm              = iniparser_getint(pardict,"CLM:nprocs",0);
  clmupdate_swc         = iniparser_getint(pardict,"CLM:update_swc",1);
  clmupdate_T           = iniparser_getint(pardict,"CLM:update_T",0);
  clmupdate_texture     = iniparser_getint(pardict,"CLM:update_texture",0);
  clmprint_swc          = iniparser_getint(pardict,"CLM:print_swc",0);
  clmprint_et           = iniparser_getint(pardict,"CLM:print_et",0);
 
  /* get settings for COSMO */
  nproccosmo      = iniparser_getint(pardict,"COSMO:nprocs",0);
  dtmult_cosmo    = iniparser_getint(pardict,"COSMO:dtmult",0);

  /* get settings for data assimilation */
  string                = iniparser_getstring(pardict,"DA:outdir","");
  strcpy(outdir,string);
  nreal                 = iniparser_getint(pardict,"DA:nreal",0);
  startreal             = iniparser_getint(pardict,"DA:startreal",0);
  da_interval           = iniparser_getdouble(pardict,"DA:da_interval",1);
  stat_dumpoffset       = iniparser_getint(pardict,"DA:stat_dumpoffset",0);
  screen_wrapper        = iniparser_getint(pardict,"DA:screen_wrapper",1);
  point_obs             = iniparser_getint(pardict,"DA:point_obs",1);
  obs_interp_switch     = iniparser_getint(pardict,"DA:obs_interp_switch",0);
  total_steps = (int) (t_sim/da_interval);
  tstartcycle = (int) (t_start/da_interval);

  /* print inputs / debug output for data assimilation settings */
  if (mype_world == 0) {
    if (screen_wrapper > 0) {
      printf("TSMP-PDAF-WRAPPER read_enkfpar: [DA]\n");
      printf("TSMP-PDAF-WRAPPER ------------------\n");
      printf("TSMP-PDAF-WRAPPER t_sim = %lf | da_interval = %lf | total_steps = %d\n",t_sim,da_interval,total_steps);
      printf("TSMP-PDAF-WRAPPER nreal = %d | n_modeltasks = %d\n",nreal,n_modeltasks);
    }
  }

  /* Check: `nreal` must be equal to n_modeltasks */
  if (nreal != n_modeltasks) {
    printf("Error: nreal must be equal to n_modeltasks.\n");
    exit(1);
  }

  /* Check: `point_obs` must be equal to either 0 or 1 */
  if (point_obs != 0 && point_obs != 1){
    printf("point_obs=%d\n", point_obs);
    printf("Error: point_obs must be equal to either 0 or 1.\n");
    exit(1);
  }

  /* Check: `npes_model = nprocpf + nprocclm + npproccosmo */
  if (nprocpf + nprocclm + nproccosmo != npes_model){
    printf("nprocpf=%d\n", nprocpf);
    printf("nprocclm=%d\n", nprocclm);
    printf("nproccosmo=%d\n", nproccosmo);
    printf("npes_model=%d\n", npes_model);
    printf("Error:  nprocpf + nprocclm + npproccosmo must be equal to npes_model.\n");
    exit(1);
  }

  /* Assign model specifier (0=clm, 1=parflow, 2=cosmo) */
  /* Order: ParFlow, CLM, COSMO */
  if (mype_model < nprocpf) {
    model = 1;
  }
  else if(mype_model < (nprocclm+nprocpf)){
    model = 0;
  }
  else{
    model = 2;
  }

  /* MPI: Get size and rank in COMM_WORLD */
  /* define number of first model realisation (for input/output filenames) */
  coupcol = task_id - 1 + startreal;
  if (screen_wrapper > 1) {
    printf("TSMP-PDAF-WRAPPER mype(w)=%d: coupcol, task_id = %d, %d\n", mype_world, coupcol,task_id);
  }

  /* create instance specific input file for ParFLow and CLM*/
  //sprintf(pfinfile ,"%s_%05d",pfinfile,coupcol);
  if((strlen(pfprefixin)) == 0){
    sprintf(pfinfile,"%s_%05d",pfproblemname,coupcol);
  }else{
    sprintf(pfinfile,"%s_%05d/%s_%05d",pfprefixin,coupcol,pfproblemname,coupcol);
  }
  sprintf(clminfile,"%s_%05d",clminfile,coupcol);
  oasprefixno  = coupcol;
  clmprefixlen = (int)strlen(clminfile);

  /* create output filenames for ParFlow */
  if((strlen(outdir)) == 0){
    strcpy(pfoutfile_stat,pfproblemname);
    if((strlen(pfprefixout))==0){
      sprintf(pfoutfile_ens,"%s_%05d",pfproblemname,coupcol);
    }else{
      sprintf(pfoutfile_ens,"%s_%05d/%s_%05d",pfprefixout,coupcol,pfproblemname,coupcol);
    }
  }else{
    sprintf(pfoutfile_stat,"%s/%s",outdir,pfproblemname);
    if((strlen(pfprefixout))==0){
      sprintf(pfoutfile_ens,"%s/%s_%05d",outdir,pfproblemname,coupcol);
    }else{
      sprintf(pfoutfile_ens,"%s/%s_%05d/%s_%05d",outdir,pfprefixout,coupcol,pfproblemname,coupcol);
    }
  }

}
