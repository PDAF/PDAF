

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
read_enkfpar.c: Function for reading controle file of TerrSysMP-PDAF
-------------------------------------------------------------------------------------------*/

#include "enkf.h"
#include "iniparser.h"

void read_enkfpar(char *parname)
{
  char *string;
  dictionary *pardict;
 
  /* initialize dictionary */
  pardict = iniparser_load(parname);
 
  /* get settings for ParFlow */
  string          = iniparser_getstring(pardict,"PF:problemname", "");
  strcat(pfinfile,string);
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
  
 
  /* get settings for CLM */
  string            = iniparser_getstring(pardict,"CLM:problemname", "");
  strcat(clminfile,string);
  nprocclm          = iniparser_getint(pardict,"CLM:nprocs",0);
  clmupdate_swc     = iniparser_getint(pardict,"CLM:update_swc",1);
  clmupdate_texture = iniparser_getint(pardict,"CLM:update_texture",0);
  clmprint_swc      = iniparser_getint(pardict,"CLM:print_swc",0);
 
  /* get settings for data assimilation */
  string          = iniparser_getstring(pardict,"DA:outdir","");
  strcat(outdir,string);
  nreal           = iniparser_getint(pardict,"DA:nreal",0);
  //stat_dumpint    = iniparser_getint(pardict,"DA:stat_dumpinterval",1);
  da_interval     = iniparser_getdouble(pardict,"DA:da_interval",1);
  stat_dumpoffset = iniparser_getint(pardict,"DA:stat_dumpoffset",0);
 
  nsteps = (int) (t_end/da_interval); 
  printf("t_end = %lf | da_interval = %lf | nsteps = %d\n",t_end,da_interval,nsteps);
  

  /* get settings for COSMO */
  nproccosmo      = iniparser_getint(pardict,"COSMO:nprocs",0);
  dtmult_cosmo    = iniparser_getint(pardict,"COSMO:dtmult",0);


/* ERROR CHECKING */
  /* check if stat_dumpinterval and da_interval are synchronized */
  //if(stat_dumpint%da_interval){
  //    printf("DA:stat_dumpint should be a multiple of DA:da_interval\nPlease check your input!\nAborting...");
  //    exit(1);
  //}

  //printf("ParFlow update flag: %d\n",pf_updateflag);
  //printf("ParFlow parameter update flag: %d\n",pf_paramupdate);
}


