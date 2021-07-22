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
wrapper_tsmp.c: Wrapper functions for TerrSysMP
-------------------------------------------------------------------------------------------*/

#define GLOBAL
#include "enkf.h"
#if defined COUP_OAS_PFL || defined PARFLOW_STAND_ALONE
#include "enkf_parflow.h"
#endif
#include "wrapper_tsmp.h"

void initialize_tsmp() {
  int rank,size;
  int subrank,subsize;
  int coupcol,coupkey;
  int interpfcol,interpfkey;
  int tcycle;
  int *pfrank;
  int i,j;
  int argc = 0; char ** argv ;


  /* read parameter file for data assimilation 'enkfpf.par' */
  read_enkfpar("enkfpf.par");


  /* assign model number (0=clm, 1=parflow, 2=cosmo) */
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  coupcol = rank / (size/nreal);
  subrank = mype_model;

  /* define number of first model realisation (for input/output filenames) */
  coupcol = coupcol + startreal;

  /* CLM, ParFlow, COSMO */
  if (subrank < nprocclm) {
    model = 0;
  }
  else if(subrank < (nprocclm+nprocpf)){
    model = 1;
  }
  else{
    model = 2;
  }
  /* ParFlow, CLM, COSMO */
  if (subrank < nprocpf) {
    model = 1;
  }
  else if(subrank < (nprocclm+nprocpf)){
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


  /* initialize clm and parflow instances */
  if(model == 0) {
#if defined COUP_OAS_PFL || defined CLMSA || defined COUP_OAS_COS
    clm_init(clminfile);
#endif
  }
  if(model == 1) {
#if defined COUP_OAS_PFL || defined PARFLOW_STAND_ALONE
    enkfparflowinit(argc,argv,pfinfile);
    idx_map_subvec2state   = (int *)   malloc(enkf_subvecsize * sizeof(int));
    parflow_oasis_init(t_start,(double)da_interval);
    
    pf_statevecsize = enkf_subvecsize;
    if(pf_updateflag == 3 || pf_updateflag == 2) pf_statevecsize = pf_statevecsize * 2;
    
    pf_paramvecsize = enkf_subvecsize;
    if(pf_paramupdate == 2) pf_paramvecsize = nx_local*ny_local;
    if(pf_paramupdate == 4 || pf_paramupdate == 5) pf_paramvecsize = 2*enkf_subvecsize;
    if(pf_paramupdate == 6 || pf_paramupdate == 7) pf_paramvecsize = 3*enkf_subvecsize;
    if(pf_paramupdate == 8) pf_paramvecsize = 4*enkf_subvecsize;
    if(pf_paramupdate > 0) pf_statevecsize += pf_paramvecsize;

    subvec_p               = (double*) calloc(enkf_subvecsize,sizeof(double));
    subvec_sat             = (double*) calloc(enkf_subvecsize,sizeof(double));
    subvec_porosity        = (double*) calloc(enkf_subvecsize,sizeof(double));
    subvec_param           = (double*) calloc(pf_paramvecsize,sizeof(double));
    subvec_mean            = (double*) calloc(enkf_subvecsize,sizeof(double));
    subvec_sd              = (double*) calloc(enkf_subvecsize,sizeof(double));
    if(pf_gwmasking > 0){
    subvec_gwind           = (double*) calloc(enkf_subvecsize,sizeof(double));
    }
    if(pf_paramupdate == 4){
        dat_alpha          = (double*) calloc(enkf_subvecsize,sizeof(double));
        dat_n              = (double*) calloc(enkf_subvecsize,sizeof(double));
    }
    else if(pf_paramupdate == 5){
        dat_ksat    = (double*) calloc(enkf_subvecsize,sizeof(double));
        dat_poro    = (double*) calloc(enkf_subvecsize,sizeof(double));
    }
    else if(pf_paramupdate == 6){
        dat_ksat     = (double*) calloc(enkf_subvecsize,sizeof(double));
        dat_alpha    = (double*) calloc(enkf_subvecsize,sizeof(double));
        dat_n        = (double*) calloc(enkf_subvecsize,sizeof(double));
    }
    else if(pf_paramupdate == 7){
        dat_poro     = (double*) calloc(enkf_subvecsize,sizeof(double));
        dat_alpha    = (double*) calloc(enkf_subvecsize,sizeof(double));
        dat_n        = (double*) calloc(enkf_subvecsize,sizeof(double));
    }
    else if(pf_paramupdate == 8){
        dat_ksat     = (double*) calloc(enkf_subvecsize,sizeof(double));
        dat_poro     = (double*) calloc(enkf_subvecsize,sizeof(double));
        dat_alpha    = (double*) calloc(enkf_subvecsize,sizeof(double));
        dat_n        = (double*) calloc(enkf_subvecsize,sizeof(double));
    }
    
    pf_statevec            = (double*) calloc(pf_statevecsize,sizeof(double));
#endif
  }

  if(model == 2){
#if defined COUP_OAS_COS
    cosmo_init();
#endif
  }
}


void finalize_tsmp() {

  if(model == 0) {
#if defined COUP_OAS_PFL || defined CLMSA || defined COUP_OAS_COS
    clm_finalize();
#endif
  }

  if(model == 1) {
#if defined COUP_OAS_PFL || defined PARFLOW_STAND_ALONE
    free(subvec_p);
    free(subvec_sat);
    free(subvec_porosity);
    free(subvec_param);
    free(pf_statevec);
    enkfparflowfinalize();
#endif
  }

  if(model == 2){
#if defined COUP_OAS_COS
    cosmo_finalize();
#endif
  }
}

void integrate_tsmp() {

  if(model == 0){
#if defined COUP_OAS_PFL || defined CLMSA || defined COUP_OAS_COS
    int tsclm;
    tsclm = (int) ( (double) da_interval / dt );
    //printf("CLM: advancing (%d clm time steps)\n",tsclm);
    clm_advance(&tsclm);
    //printf("CLM: advancing finished\n",tsclm);
#endif
  }

  if(model == 1){
#if defined COUP_OAS_PFL || defined PARFLOW_STAND_ALONE
    //printf("Parflow: advancing (from %lf to %lf)\n",t_start,t_start+(double)da_interval);
    enkfparflowadvance(t_start,(double)da_interval);
   //// printf("Parflow: advancing finished\n");
//    enkf_printstatistics_pfb(pf_statevec,"state",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_ens,3);

    if(pf_printstat==1){
      MPI_Comm comm_couple_c = MPI_Comm_f2c(comm_couple);
      enkf_ensemblestatistics(pf_statevec,subvec_mean,subvec_sd,enkf_subvecsize,comm_couple_c);
      if(task_id==1 && pf_updateflag==1){
        enkf_printstatistics_pfb(subvec_mean,"press.mean",(int) (t_start/da_interval + 1 + stat_dumpoffset),pfoutfile_stat,3);
        enkf_printstatistics_pfb(subvec_sd,"press.sd",(int) (t_start/da_interval + 1 + stat_dumpoffset),pfoutfile_stat,3);
      }
      if(task_id==1 && (pf_updateflag==3 || pf_updateflag==2)){
        enkf_printstatistics_pfb(subvec_mean,"swc.mean",(int) (t_start/da_interval + 1 + stat_dumpoffset ),pfoutfile_stat,3);
        enkf_printstatistics_pfb(subvec_sd,"swc.sd",(int) (t_start/da_interval + 1 + stat_dumpoffset),pfoutfile_stat,3);
      }
    }
#endif
  }

  if(model == 2){
#if defined COUP_OAS_COS
    int tscos;
    tscos = (int) ((double) da_interval / dt);
    tscos = tscos * dtmult_cosmo;
    //printf("tscos is %d",tscos);
    cosmo_advance(&tscos);
#endif
  }

  //print_memusage((int) t_start);
  t_start += (double)da_interval;
}

#if (defined COUP_OAS_PFL || defined PARFLOW_STAND_ALONE)
void print_update_pfb(){
  if(model == 1){
    enkf_printstatistics_pfb(subvec_p,"update",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_ens,3);
  }
}
#endif



void update_tsmp(){

  /* print analysis and update parflow */
#if (defined COUP_OAS_PFL || defined PARFLOW_STAND_ALONE)
  if(model == 1){
    int i;
    double *dat;
    int do_pupd=0;

    /* check if frequency of parameter update is reached */
    do_pupd = (int) t_start/da_interval;
    do_pupd = do_pupd % pf_freq_paramupdate;
    do_pupd = !do_pupd;

    /* print updated ensemble */
    if(pf_updateflag == 3){
      dat = &pf_statevec[enkf_subvecsize];
    }else{
      dat = pf_statevec;
    }
    if(pf_printensemble == 1) enkf_printstatistics_pfb(dat,"update",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_ens,3);

    /* update Ksat */
    if(pf_paramupdate == 1 && do_pupd){
      dat = &pf_statevec[pf_statevecsize-pf_paramvecsize];

      /* dampening */
      for(i=0;i<pf_paramvecsize;i++) dat[i] = log(subvec_param[i]) + pf_dampfac_param * (dat[i] - log(subvec_param[i]));

      /* print ensemble statistics */
      if(pf_paramprintstat){
        MPI_Comm comm_couple_c = MPI_Comm_f2c(comm_couple);
        enkf_ensemblestatistics(dat,subvec_mean,subvec_sd,pf_paramvecsize,comm_couple_c);
        if(task_id==1){
          enkf_printstatistics_pfb(subvec_mean,"param.ksat.mean",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
          enkf_printstatistics_pfb(subvec_sd,"param.ksat.sd",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
        }
      }

      /* backtransform updated K values */
      for(i=0;i<pf_paramvecsize;i++) dat[i] = exp(dat[i]);

      /* print updated K values */
      if(pf_paramprintensemble) enkf_printstatistics_pfb(dat,"update.param.ksat",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_ens,3);
    }


    /* update Mannings */
    if(pf_paramupdate == 2 && do_pupd){
      dat = &pf_statevec[pf_statevecsize-pf_paramvecsize];

      /* dampening */
      for(i=0;i<pf_paramvecsize;i++) dat[i] = log10(subvec_param[i]) + pf_dampfac_param * (dat[i] - log10(subvec_param[i]));

      /* print ensemble statistics */
      if(pf_paramprintstat){
        MPI_Comm comm_couple_c = MPI_Comm_f2c(comm_couple);
        enkf_ensemblestatistics(dat,subvec_mean,subvec_sd,pf_paramvecsize,comm_couple_c);
        if(task_id==1){
          enkf_printstatistics_pfb(subvec_mean,"param.mannings.mean",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,2);
          enkf_printstatistics_pfb(subvec_sd,"param.mannings.sd",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,2);
        }
      }

      /* backtransform updated mannings values */
      dat = &pf_statevec[pf_statevecsize-pf_paramvecsize];
      for(i=0;i<pf_paramvecsize;i++) dat[i] = pow(10,dat[i]);

      /* print updated mannings values */
      if(pf_paramprintensemble){
        char fprefix [200];
        char fsuffix [10];
        //sprintf(fprefix,"%s/%s.%s",outdir,pfinfile,"update.mannings");
        sprintf(fprefix,"%s.%s",pfoutfile_ens,"update.mannings");
        sprintf(fsuffix,"%05d",(int) (t_start/da_interval + stat_dumpoffset));
        enkf_printmannings(fprefix,fsuffix);
      }
    }

    /* update porosity */
    if(pf_paramupdate == 3 && do_pupd){
      dat = &pf_statevec[pf_statevecsize-pf_paramvecsize];

      /* dampening */
      for(i=0;i<pf_paramvecsize;i++) dat[i] = subvec_param[i] + pf_dampfac_param * (dat[i] - subvec_param[i]);

      /* print ensemble statistics */
      if(pf_paramprintstat){
        MPI_Comm comm_couple_c = MPI_Comm_f2c(comm_couple);
        enkf_ensemblestatistics(dat,subvec_mean,subvec_sd,pf_paramvecsize,comm_couple_c);
        if(task_id == 1){
          enkf_printstatistics_pfb(subvec_mean,"param.poro.mean", (int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
          enkf_printstatistics_pfb(subvec_sd,"param.poro.sd", (int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
        }
      }

      /* print updated porosity values */
      if(pf_paramprintensemble) enkf_printstatistics_pfb(dat, "update.param.poro", (int) (t_start/da_interval + stat_dumpoffset),pfoutfile_ens,3);

    }
    
    /* update van Genuchten */
    if(pf_paramupdate == 4 && do_pupd){
        dat = &pf_statevec[pf_statevecsize-pf_paramvecsize];

        /* dampening */
        int alpha_counter = 0;
        int n_counter = 0;
        for(i=0;i<pf_paramvecsize;i++){
          if((i%2)==0){
            dat[i] = log(subvec_param[i]) + pf_dampfac_param * (dat[i] - log(subvec_param[i]));
            dat_alpha[alpha_counter] = dat[i];
            alpha_counter++;
          }else{
            dat[i] = subvec_param[i] + pf_dampfac_param * (dat[i] - subvec_param[i]);
            dat_n[n_counter] = dat[i];
            n_counter++;      
           }
        }

        /* print ensemble statistics */
        if(pf_paramprintstat){
          MPI_Comm comm_couple_c = MPI_Comm_f2c(comm_couple);
          enkf_ensemblestatistics(dat_alpha,subvec_mean,subvec_sd,pf_paramvecsize/2,comm_couple_c);
          if(task_id==1){
            enkf_printstatistics_pfb(subvec_mean,"param.alpha.mean",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
            enkf_printstatistics_pfb(subvec_sd,"param.alpha.sd",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
          }
          enkf_ensemblestatistics(dat_n,subvec_mean,subvec_sd,pf_paramvecsize/2,comm_couple_c);
          if(task_id==1){
            enkf_printstatistics_pfb(subvec_mean,"param.n.mean",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
            enkf_printstatistics_pfb(subvec_sd,"param.n.sd",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
          }
        }

        /* backtransform updated alpha values */
        alpha_counter = 0;
        for(i=0;i<pf_paramvecsize;i++){
            if((i%2)==0){
                dat[i] = exp(dat[i]);
                dat_alpha[alpha_counter] = dat[i];
                alpha_counter++;
            }
        }

        /* print updated van Genuchten values */
        if(pf_paramprintensemble){
            enkf_printstatistics_pfb(dat_alpha,"update.param.alpha",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_ens,3);
            enkf_printstatistics_pfb(dat_n,"update.param.n",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_ens,3);
        }
    }
    
    /* update hydraulic conductivity and porosity */
    if(pf_paramupdate == 5 && do_pupd){
        dat = &pf_statevec[pf_statevecsize-pf_paramvecsize];
        
        /* dampening */
        int ksat_counter = 0;
        int poro_counter = 0;
        for(i=0;i<pf_paramvecsize;i++){
            if((i%2)==0){
                dat[i] = log(subvec_param[i]) + pf_dampfac_param * (dat[i] - log(subvec_param[i]));
                dat_ksat[ksat_counter] = dat[i];
                ksat_counter++;
            }else{
                dat[i] = subvec_param[i] + pf_dampfac_param * (dat[i] - subvec_param[i]);
                dat_poro[poro_counter] = dat[i];
                poro_counter++;
            }
        }
        
        /* print ensemble statistics */
        if(pf_paramprintstat){
            MPI_Comm comm_couple_c = MPI_Comm_f2c(comm_couple);
            enkf_ensemblestatistics(dat_ksat,subvec_mean,subvec_sd,pf_paramvecsize/2,comm_couple_c);
            if(task_id==1){
              enkf_printstatistics_pfb(subvec_mean,"param.ksat.mean",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
              enkf_printstatistics_pfb(subvec_sd,"param.ksat.sd",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
            }
            enkf_ensemblestatistics(dat_poro,subvec_mean,subvec_sd,pf_paramvecsize/2,comm_couple_c);
            if(task_id==1){
              enkf_printstatistics_pfb(subvec_mean,"param.poro.mean",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
              enkf_printstatistics_pfb(subvec_sd,"param.poro.sd",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
            }
        }
        
        /* backtransform updated ksat values */
        ksat_counter = 0;
        for(i=0;i<pf_paramvecsize;i++){
            if((i%2)==0){
                dat[i] = exp(dat[i]);
                dat_ksat[ksat_counter] = dat[i];
                ksat_counter++;
            }
        }

        /* print updated parameter values */
        if(pf_paramprintensemble){ 
            enkf_printstatistics_pfb(dat_ksat,"update.param.ksat",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_ens,3);
            enkf_printstatistics_pfb(dat_poro,"update.param.poro",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_ens,3);
        }
    }
    
    /* update hydraulic conductivity and van Genuchten parameters */
    if(pf_paramupdate == 6 && do_pupd){
        dat = &pf_statevec[pf_statevecsize-pf_paramvecsize];
        
        /* dampening */
        int ksat_counter = 0;
        int alpha_counter = 0;
        int n_counter = 0;
        for(i=0;i<pf_paramvecsize;i=i+3){
            dat[i] = log(subvec_param[i]) + pf_dampfac_param * (dat[i] - log(subvec_param[i]));
            dat_ksat[ksat_counter] = dat[i];
            ksat_counter++;
            dat[i+1] = log(subvec_param[i+1]) + pf_dampfac_param * (dat[i+1] - log(subvec_param[i+1]));
            dat_alpha[alpha_counter] = dat[i+1];
            alpha_counter++;
            dat[i+2] = subvec_param[i+2] + pf_dampfac_param * (dat[i+2] - subvec_param[i+2]);
            dat_n[n_counter] = dat[i+2];
            n_counter++;
        }
        
        /* print ensemble statistics */
        if(pf_paramprintstat){
            MPI_Comm comm_couple_c = MPI_Comm_f2c(comm_couple);
            enkf_ensemblestatistics(dat_ksat,subvec_mean,subvec_sd,pf_paramvecsize/3,comm_couple_c);
            if(task_id==1){
              enkf_printstatistics_pfb(subvec_mean,"param.ksat.mean",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
              enkf_printstatistics_pfb(subvec_sd,"param.ksat.sd",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
            }
            enkf_ensemblestatistics(dat_alpha,subvec_mean,subvec_sd,pf_paramvecsize/3,comm_couple_c);
            if(task_id==1){
              enkf_printstatistics_pfb(subvec_mean,"param.alpha.mean",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
              enkf_printstatistics_pfb(subvec_sd,"param.alpha.sd",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
            }
            enkf_ensemblestatistics(dat_n,subvec_mean,subvec_sd,pf_paramvecsize/3,comm_couple_c);
            if(task_id==1){
              enkf_printstatistics_pfb(subvec_mean,"param.n.mean",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
              enkf_printstatistics_pfb(subvec_sd,"param.n.sd",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
            }
        }
        
        /* backtransform updated ksat values */
        ksat_counter = 0;
        alpha_counter = 0;
        for(i=0;i<pf_paramvecsize;i=i+3){
            dat[i] = exp(dat[i]);
            dat_ksat[ksat_counter] = dat[i];
            ksat_counter++;
            dat[i+1] = exp(dat[i+1]);
            dat_alpha[alpha_counter] = dat[i+1];
            alpha_counter++;
        }

        /* print updated parameter values */
        if(pf_paramprintensemble){ 
            enkf_printstatistics_pfb(dat_ksat,"update.param.ksat",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_ens,3);
            enkf_printstatistics_pfb(dat_alpha,"update.param.alpha",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_ens,3);
            enkf_printstatistics_pfb(dat_n,"update.param.n",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_ens,3);
        }
    }
    
    /* update porosity and van Genuchten parameters */
    if(pf_paramupdate == 7 && do_pupd){
        dat = &pf_statevec[pf_statevecsize-pf_paramvecsize];
        
        /* dampening */
        int poro_counter = 0;
        int alpha_counter = 0;
        int n_counter = 0;
        for(i=0;i<pf_paramvecsize;i=i+3){
            dat[i] = subvec_param[i] + pf_dampfac_param * (dat[i] - subvec_param[i]);
            dat_poro[poro_counter] = dat[i];
            poro_counter++;
            dat[i+1] = log(subvec_param[i+1]) + pf_dampfac_param * (dat[i+1] - log(subvec_param[i+1]));
            dat_alpha[alpha_counter] = dat[i+1];
            alpha_counter++;
            dat[i+2] = subvec_param[i+2] + pf_dampfac_param * (dat[i+2] - subvec_param[i+2]);
            dat_n[n_counter] = dat[i+2];
            n_counter++;
        }
        
        /* print ensemble statistics */
        if(pf_paramprintstat){
            MPI_Comm comm_couple_c = MPI_Comm_f2c(comm_couple);
            enkf_ensemblestatistics(dat_poro,subvec_mean,subvec_sd,pf_paramvecsize/3,comm_couple_c);
            if(task_id==1){
              enkf_printstatistics_pfb(subvec_mean,"param.poro.mean",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
              enkf_printstatistics_pfb(subvec_sd,"param.poro.sd",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
            }
            enkf_ensemblestatistics(dat_alpha,subvec_mean,subvec_sd,pf_paramvecsize/3,comm_couple_c);
            if(task_id==1){
              enkf_printstatistics_pfb(subvec_mean,"param.alpha.mean",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
              enkf_printstatistics_pfb(subvec_sd,"param.alpha.sd",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
            }
            enkf_ensemblestatistics(dat_n,subvec_mean,subvec_sd,pf_paramvecsize/3,comm_couple_c);
            if(task_id==1){
              enkf_printstatistics_pfb(subvec_mean,"param.n.mean",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
              enkf_printstatistics_pfb(subvec_sd,"param.n.sd",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
            }
        }
        
        /* backtransform updated ksat values */
        alpha_counter = 0;
        for(i=1;i<pf_paramvecsize;i=i+3){
            dat[i] = exp(dat[i]);
            dat_alpha[alpha_counter] = dat[i];
            alpha_counter++;
        }

        /* print updated parameter values */
        if(pf_paramprintensemble){ 
            enkf_printstatistics_pfb(dat_poro,"update.param.poro",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_ens,3);
            enkf_printstatistics_pfb(dat_alpha,"update.param.alpha",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_ens,3);
            enkf_printstatistics_pfb(dat_n,"update.param.n",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_ens,3);
        }
    }
    
    /* update hydraulic conductivity, porosity and van Genuchten parameters */
    if(pf_paramupdate == 8 && do_pupd){
        dat = &pf_statevec[pf_statevecsize-pf_paramvecsize];
        
        /* dampening */
        int ksat_counter = 0;
        int poro_counter = 0;
        int alpha_counter = 0;
        int n_counter = 0;
        for(i=0;i<pf_paramvecsize;i=i+4){
            dat[i] = log(subvec_param[i]) + pf_dampfac_param * (dat[i] - log(subvec_param[i]));
            dat_ksat[ksat_counter] = dat[i];
            ksat_counter++;
            dat[i+1] = subvec_param[i+1] + pf_dampfac_param * (dat[i+1] - subvec_param[i+1]);
            dat_poro[poro_counter] = dat[i+1];
            poro_counter++;
            dat[i+2] = log(subvec_param[i+2]) + pf_dampfac_param * (dat[i+2] - log(subvec_param[i+2]));
            dat_alpha[alpha_counter] = dat[i+2];
            alpha_counter++;
            dat[i+3] = subvec_param[i+3] + pf_dampfac_param * (dat[i+3] - subvec_param[i+3]);
            dat_n[n_counter] = dat[i+3];
            n_counter++;
        }
        
        /* print ensemble statistics */
        if(pf_paramprintstat){
            MPI_Comm comm_couple_c = MPI_Comm_f2c(comm_couple);
            enkf_ensemblestatistics(dat_ksat,subvec_mean,subvec_sd,pf_paramvecsize/4,comm_couple_c);
            if(task_id==1){
              enkf_printstatistics_pfb(subvec_mean,"param.ksat.mean",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
              enkf_printstatistics_pfb(subvec_sd,"param.ksat.sd",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
            }
            enkf_ensemblestatistics(dat_poro,subvec_mean,subvec_sd,pf_paramvecsize/4,comm_couple_c);
            if(task_id==1){
              enkf_printstatistics_pfb(subvec_mean,"param.poro.mean",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
              enkf_printstatistics_pfb(subvec_sd,"param.poro.sd",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
            }
            enkf_ensemblestatistics(dat_alpha,subvec_mean,subvec_sd,pf_paramvecsize/4,comm_couple_c);
            if(task_id==1){
              enkf_printstatistics_pfb(subvec_mean,"param.alpha.mean",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
              enkf_printstatistics_pfb(subvec_sd,"param.alpha.sd",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
            }
            enkf_ensemblestatistics(dat_n,subvec_mean,subvec_sd,pf_paramvecsize/4,comm_couple_c);
            if(task_id==1){
              enkf_printstatistics_pfb(subvec_mean,"param.n.mean",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
              enkf_printstatistics_pfb(subvec_sd,"param.n.sd",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_stat,3);
            }
        }
        
        /* backtransform updated ksat values */
        ksat_counter = 0;
        alpha_counter = 0;
        for(i=0;i<pf_paramvecsize;i=i+4){
            dat[i] = exp(dat[i]);
            dat_ksat[ksat_counter] = dat[i];
            ksat_counter++;
            dat[i+2] = exp(dat[i+2]);
            dat_alpha[alpha_counter] = dat[i+2];
            alpha_counter++;
        }

        /* print updated parameter values */
        if(pf_paramprintensemble){ 
            enkf_printstatistics_pfb(dat_ksat,"update.param.ksat",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_ens,3);
            enkf_printstatistics_pfb(dat_poro,"update.param.poro",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_ens,3);
            enkf_printstatistics_pfb(dat_alpha,"update.param.alpha",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_ens,3);
            enkf_printstatistics_pfb(dat_n,"update.param.n",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_ens,3);
        }
    }
    
    update_parflow(do_pupd);
  }
#endif

}

