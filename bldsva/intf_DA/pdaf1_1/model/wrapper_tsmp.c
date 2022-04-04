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

/*-------------------------------------------------------------------------*/
/**
  @author   Wolfgang Kurtz, Guowei He
  @brief    Initialization of component models for for TSMP-PDAF.

  1. read parameter file for data assimilation 'enkfpf.par'
  2. initialize clm, parflow and cosmo instances
 */
/*--------------------------------------------------------------------------*/
void initialize_tsmp() {
  int argc = 0; char ** argv ;	/* Dummy command line arguments for amps */


  /* read parameter file for data assimilation 'enkfpf.par' */
  read_enkfpar("enkfpf.par");

  /* initialize clm, parflow and cosmo instances */
  if(model == 0) {
#if defined COUP_OAS_PFL || defined CLMSA || defined COUP_OAS_COS
    /* enkf_clm.F90 */
    clm_init(clminfile);
#endif
  }
  if(model == 1) {
#if defined COUP_OAS_PFL || defined PARFLOW_STAND_ALONE
    /* enkf_parflow.c */
    enkfparflowinit(argc,argv,pfinfile);
    parflow_oasis_init(t_start,(double)da_interval);
#endif
  }

  if(model == 2){
#if defined COUP_OAS_COS
    /* enkf_cosmo.F90 */
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
    free(subvec_mean);
    free(subvec_sd);
    free(subvec_param_mean);
    free(subvec_param_sd);
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

/*-------------------------------------------------------------------------*/
/**
  @author   Wolfgang Kurtz, Guowei He
  @brief    Integration/Simulation of component models.

  Depending on value of `model` for the current PE:
  1. Integration of CLM
  or
  1. Integration of ParFlow
  or
  1. Integration of COSMO

  - Debug output before and after integrating a given component model
    is written when `screen_wrapper > 1`.
  - For ParFlow, ensemble statistics are computed and printed to PFB files
    according to input `PF:printstat`
  - For CLM and COSMO, the number of time steps is computed before integrating.
  - For COSMO a multiplier is applied from input `COSMO:dtmult`.

  2. Advance `t_start` by `da_interval`.
 */
/*--------------------------------------------------------------------------*/
void integrate_tsmp() {

  /* CLM */
  if(model == 0){
#if defined COUP_OAS_PFL || defined CLMSA || defined COUP_OAS_COS
    /* Number of time steps for CLM */
    int tsclm;
    tsclm = (int) ( (double) da_interval / dt );

    /* Debug output */
    if (screen_wrapper > 1) {
      printf("TSMP-PDAF-WRAPPER mype(w)=%d: CLM: advancing (%d clm time steps)\n",mype_world, tsclm);
    }

    /* Integrate CLM */
    clm_advance(&tsclm);

    /* Debug output */
    if (screen_wrapper > 1) {
      printf("TSMP-PDAF-WRAPPER mype(w)=%d: CLM: advancing finished\n", mype_world);
    }

#endif
  }

  /* ParFlow */
  if(model == 1){
#if defined COUP_OAS_PFL || defined PARFLOW_STAND_ALONE
    /* Debug output */
    if (screen_wrapper > 1) {
      printf("TSMP-PDAF-WRAPPER mype(w)=%d: Parflow: advancing (from %lf to %lf)\n",mype_world,t_start,t_start+(double)da_interval);
    }

    /* Integrate ParFlow */
    enkfparflowadvance(tcycle, t_start,(double)da_interval);

    /* Debug output */
    if (screen_wrapper > 1) {
      printf("TSMP-PDAF-WRAPPER mype(w)=%d: Parflow: advancing finished\n", mype_world);
    }

    /* Print ensemble statistics to PFB */
    if(pf_printstat==1){
      printstat_parflow();
    }
#endif
  }

  /* COSMO */
  if(model == 2){
#if defined COUP_OAS_COS

    /* Number of time steps for COSMO */
    int tscos;
    tscos = (int) ((double) da_interval / dt);
    tscos = tscos * dtmult_cosmo; /* Multiplier read from input */

    /* Debug output */
    if (screen_wrapper > 1) {
      printf("TSMP-PDAF-WRAPPER mype(w)=%d: COSMO: tscos is %d",tscos);
    }

    /* Integrate COSMO */
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

#if defined CLMSA
  if((model == tag_model_clm) && (clmupdate_swc != 0)){
    update_clm();
    print_update_clm(tcycle, total_steps);
  }
#endif

  /* print analysis and update parflow */
#if (defined COUP_OAS_PFL || defined PARFLOW_STAND_ALONE)
  if(model == 1){
    int i;
    double *dat;
    int do_pupd=0;

    /* print updated ensemble */
    if(pf_updateflag == 3){
      dat = &pf_statevec[enkf_subvecsize];
    }else{
      dat = pf_statevec;
    }
    if(pf_printensemble == 1) enkf_printstatistics_pfb(dat,"update",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_ens,3);


    /* check if frequency of parameter update is reached */
    do_pupd = (int) t_start/da_interval;
    do_pupd = do_pupd % pf_freq_paramupdate;
    do_pupd = !do_pupd;

    /* update Ksat */
    if(pf_paramupdate == 1 && do_pupd){
      dat = &pf_statevec[pf_statevecsize-pf_paramvecsize];

      /* dampening */
      for(i=0;i<pf_paramvecsize;i++) dat[i] = log10(subvec_param[i]) + pf_dampfac_param * (dat[i] - log10(subvec_param[i]));

      /* print ensemble statistics */
      if(pf_paramprintstat){
	printstat_param_parflow(dat, "param", 3);
      }

      /* backtransform updated K values */
      for(i=0;i<pf_paramvecsize;i++) dat[i] = pow(10,dat[i]);

      /* print updated K values */
      if(pf_paramprintensemble) enkf_printstatistics_pfb(dat,"update.param",(int) (t_start/da_interval + stat_dumpoffset),pfoutfile_ens,3);
    }


    /* update Mannings */
    if(pf_paramupdate == 2 && do_pupd){
      dat = &pf_statevec[pf_statevecsize-pf_paramvecsize];

      /* dampening */
      for(i=0;i<pf_paramvecsize;i++) dat[i] = log10(subvec_param[i]) + pf_dampfac_param * (dat[i] - log10(subvec_param[i]));

      /* print ensemble statistics */
      if(pf_paramprintstat){
	printstat_param_parflow(dat, "param", 2);
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

    update_parflow(do_pupd);

    /* print updated mannings values */
    //if(pf_paramupdate == 2){
    //  char fprefix [200];
    //  char fsuffix [10];
    //  sprintf(fprefix,"%s/%s.%s",outdir,pfinfile,"update.mannings");
    //  sprintf(fsuffix,"%05d",(int) (t_start/da_interval + stat_dumpoffset));
    //  enkf_printmannings(fprefix,fsuffix);
    //}
  }
#endif

  // print et statistics
#if !defined PARFLOW_STAND_ALONE
  if(model == tag_model_clm && clmprint_et == 1){
    write_clm_statistics(tcycle, total_steps);
  }
#endif

  //  !print *,"Finished update_tsmp()"


}
