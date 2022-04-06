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
read_enkfpar.c: Function for reading controle file of TerrSysMP-PDAF
-------------------------------------------------------------------------------------------*/

#include "enkf.h"
#include "iniparser.h"

int countDigit(int n)
{
	if (n == 0)
		return -1;
	return 1 + countDigit(n / 10);
}

void read_enkfpar(char *parname)
{
  char *string;
  dictionary *pardict;
  int len;

  /* int rank,size; */
  /* int subrank; */
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
  t_end                 = iniparser_getdouble(pardict,"PF:endtime",0);
  dt                    = iniparser_getdouble(pardict,"PF:dt",0);
  pf_updateflag         = iniparser_getint(pardict,"PF:updateflag",1);
  pf_paramupdate        = iniparser_getint(pardict,"PF:paramupdate",0);
  pf_aniso_perm_y       = iniparser_getdouble(pardict,"PF:aniso_perm_y",1);
  pf_aniso_perm_z       = iniparser_getdouble(pardict,"PF:aniso_perm_z",1);
  pf_printensemble      = iniparser_getint(pardict,"PF:printensemble",1);
  pf_printstat          = iniparser_getint(pardict,"PF:printstat",1);
  pf_paramprintensemble = iniparser_getint(pardict,"PF:paramprintensemble",1);
  pf_paramprintstat     = iniparser_getint(pardict,"PF:paramprintstat",1);
  pf_olfmasking         = iniparser_getint(pardict,"PF:olfmasking",0);
  pf_gwmasking          = iniparser_getint(pardict,"PF:gwmasking",0);
  pf_printgwmask        = iniparser_getint(pardict,"PF:printgwmask",0);
  pf_dampfac_param      = iniparser_getdouble(pardict,"PF:dampingfactor_param",1.0);
  pf_freq_paramupdate   = iniparser_getint(pardict,"PF:paramupdate_frequency",1);
  
  /* get settings for CLM */
  string                = iniparser_getstring(pardict,"CLM:problemname", "");
  strcpy(clminfile,string);
  nprocclm              = iniparser_getint(pardict,"CLM:nprocs",0);
  clmupdate_swc         = iniparser_getint(pardict,"CLM:update_swc",1);
  clmupdate_texture     = iniparser_getint(pardict,"CLM:update_texture",0);
  clmprint_swc          = iniparser_getint(pardict,"CLM:print_swc",0);
  clmprint_et           = iniparser_getint(pardict,"CLM:print_et",0);
 
  /* get settings for data assimilation */
  string                = iniparser_getstring(pardict,"DA:outdir","");
  strcpy(outdir,string);
  nreal                 = iniparser_getint(pardict,"DA:nreal",0);
  startreal             = iniparser_getint(pardict,"DA:startreal",0);
  //stat_dumpint          = iniparser_getint(pardict,"DA:stat_dumpinterval",1);
  da_interval           = iniparser_getdouble(pardict,"DA:da_interval",1);
  stat_dumpoffset       = iniparser_getint(pardict,"DA:stat_dumpoffset",0);
  screen_wrapper        = iniparser_getint(pardict,"DA:screen_wrapper",1);
  point_obs             = iniparser_getint(pardict,"DA:point_obs",1);
  len = countDigit(point_obs);
  if (len > 1)
    point_obs=1;
  total_steps = (int) (t_end/da_interval);

  /* print inputs / debug output for data assimilation settings */
  if (mype_world == 0) {
    if (screen_wrapper > 0) {
      printf("TSMP-PDAF-WRAPPER read_enkfpar: [DA]\n");
      printf("TSMP-PDAF-WRAPPER ------------------\n");
      printf("TSMP-PDAF-WRAPPER t_end = %lf | da_interval = %lf | total_steps = %d\n",t_end,da_interval,total_steps);
      printf("TSMP-PDAF-WRAPPER nreal = %d | n_modeltasks = %d\n",nreal,n_modeltasks);
      if (nreal != n_modeltasks) {
	printf("Error: nreal must be equal to n_modeltasks.\n");
	exit(1);
      }
    }
  }

  /* get settings for COSMO */
  nproccosmo      = iniparser_getint(pardict,"COSMO:nprocs",0);
  dtmult_cosmo    = iniparser_getint(pardict,"COSMO:dtmult",0);


  /* MPI_Comm_size(MPI_COMM_WORLD,&size); */
  /* MPI_Comm_rank(MPI_COMM_WORLD,&rank); */
  /* coupcol = task_id - 1; */
  /* subrank = mype_model; */

  /* MPI: Get size and rank in COMM_WORLD */
  /* define number of first model realisation (for input/output filenames) */
  /* startreal: read from input in read_enkfpar */
  coupcol = task_id - 1 + startreal;
  if (screen_wrapper > 1) {
    printf("TSMP-PDAF-WRAPPER mype(w)=%d: coupcol, task_id = %d, %d\n", mype_world, coupcol,task_id);
    /* printf("DBG: size, npes_world = %d, %d\n",size,npes_world); */
    /* printf("DBG: rank, mype_world = %d, %d\n",rank,mype_world); */
    /* printf("DBG: mype_model, npes_model = %d, %d\n",mype_model,npes_model); */
  }

  /* CLM, ParFlow, COSMO */
  /* assign model number (0=clm, 1=parflow, 2=cosmo) */
  if (mype_model < nprocclm) {
    model = 0;
  }
  else if(mype_model < (nprocclm+nprocpf)){
    model = 1;
  }
  else{
    model = 2;
  }
  /* ParFlow, CLM, COSMO */
  if (mype_model < nprocpf) {
    model = 1;
  }
  else if(mype_model < (nprocclm+nprocpf)){
    model = 0;
  }
  else{
    model = 2;
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

  /* Set variables from input */


/* ERROR CHECKING */
  /* check if stat_dumpinterval and da_interval are synchronized */
  //if(stat_dumpint%da_interval){
  //    printf("DA:stat_dumpint should be a multiple of DA:da_interval\nPlease check your input!\nAborting...");
  //    exit(1);
  //}

  //printf("ParFlow update flag: %d\n",pf_updateflag);
  //printf("ParFlow parameter update flag: %d\n",pf_paramupdate);
}
