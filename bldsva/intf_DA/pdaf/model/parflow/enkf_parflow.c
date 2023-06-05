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
enkf_parflow.c: Wrapper functions for ParFlow
-------------------------------------------------------------------------------------------*/
#include "enkf.h"
#include "parflow.h"
#include "solver.h"

#include "enkf_parflow.h"
#include <string.h>
#include "problem_saturationtopressure.h"

//#include <spi/include/l1p/sprefetch.h>

amps_ThreadLocalDcl(PFModule *, Solver_module);
amps_ThreadLocalDcl(PFModule *, solver);
amps_ThreadLocalDcl(Vector *, evap_trans);
amps_ThreadLocalDcl(Vector *, vdummy_3d);
amps_ThreadLocalDcl(Vector *, vdummy_2d);
amps_ThreadLocalDcl(PFModule *, problem);

//ProblemData *GetProblemDataRichards(PFModule *this_module);
//Problem *GetProblemRichards(PFModule *this_module);
//PFModule *GetICPhasePressureRichards(PFModule *this_module);

Vector *GetPressureRichards(PFModule *this_module);
Vector *GetSaturationRichards(PFModule *this_module);
Vector *GetDensityRichards(PFModule *this_module);
int     GetEvapTransFile(PFModule *this_module);
char   *GetEvapTransFilename(PFModule *this_module);

void init_idx_map_subvec2state(Vector *pf_vector) {
	Grid *grid = VectorGrid(pf_vector);

	int sg;
	double *tmpdat;

	// allocate x, y z coords
	int num = enkf_subvecsize;

   // hcp param update conditional statement
   // we need to indicate the physical coordinates
   // of the parameter (K_sat) in the x/ycoord if
   // it is included in the state vector for
   // localization purposes.
	if( pf_paramupdate == 1 )
	   num *= 2;

	xcoord = (double *) malloc(num * sizeof(double));
	ycoord = (double *) malloc(num * sizeof(double));
	zcoord = (double *) malloc(num * sizeof(double));
	//tmpdat = (double *) malloc(enkf_subvecsize * sizeof(double));

	// copy dz_mult to double
	ProblemData * problem_data = GetProblemDataRichards(solver);
	Vector      * dz_mult      = ProblemDataZmult(problem_data);
	Problem     * problem      = GetProblemRichards(solver);

    //Vector *pressure_in = (Vector*) GetPressureRichards(amps_ThreadLocal(solver));
    //Vector *pressure_in = (solver->instance_xtra->pressure);
    //Vector *pressure_in = vdummy_3d;
    //Vector *pressure_in = GetPressureRichards2(amps_ThreadLocal(solver));
    //Vector *pressure_in = NULL;
    //GetPressureRichards3(amps_ThreadLocal(solver),&pressure_in);

	//Grid *grid = VectorGrid(pressure_in);

	//PFModule * dz_mult_module = ProblemdzScale(problem);
	//char * name;
	//name = "dzScale.nzListNumber";

	//int num_dz = GetDouble(name);
	//int ir;
	//double values[num_dz];
	//char key[IDB_MAX_KEY_LEN];
	//for (ir = 0; ir < num_dz; ir++) {
	//	sprintf(key, "Cell.%d.dzScale.Value", ir);
	//	values[ir] = GetDouble(key);
	//}
	//PF2ENKF(dz_mult, zcoord);

	ForSubgridI(sg, GridSubgrids(grid))
	{
		Subgrid *subgrid = GridSubgrid(grid, sg);

		int ix = SubgridIX(subgrid);
		int iy = SubgridIY(subgrid);
		int iz = SubgridIZ(subgrid);

		int nx = SubgridNX(subgrid);
		int ny = SubgridNY(subgrid);
		int nz = SubgridNZ(subgrid);

		nx_glob = BackgroundNX(GlobalsBackground);
		ny_glob = BackgroundNY(GlobalsBackground);
	       	nz_glob = BackgroundNZ(GlobalsBackground);

		int i, j, k;
		int counter = 0;

		for (k = iz; k < iz + nz; k++) {
			for (j = iy; j < iy + ny; j++) {
				for (i = ix; i < ix + nx; i++) {
					idx_map_subvec2state[counter] = nx_glob * ny_glob * k	+ nx_glob * j + i;
					idx_map_subvec2state[counter] += 1; // wolfgang's fix for C -> Fortran index
					//xcoord[counter] = i * SubgridDX(subgrid) + 0.5*SubgridDX(subgrid); //SubgridX(subgrid) ;//+ i * SubgridDX(subgrid);
					//ycoord[counter] = j * SubgridDY(subgrid) + 0.5*SubgridDY(subgrid); //SubgridY(subgrid) ;//+ j * SubgridDY(subgrid);
					xcoord[counter] =   i;
					ycoord[counter] =   j;
					zcoord[counter] =   k;
					//zcoord[counter] = SubgridZ(subgrid) + k * SubgridDZ(subgrid)*values[k];
                                        //tmpdat[counter] = (double)idx_map_subvec2state[counter];
					counter++;

				}
			}
		}

      //  hcp paramupdate
      //  here we indicate the physical coordinates of the
      //  parameters according to their addresses in the
      //  state vector.
      if( pf_paramupdate == 1 )
      {
         for( i = 0; i < enkf_subvecsize; i++ ) {
            xcoord[enkf_subvecsize + i] = xcoord[i];
            ycoord[enkf_subvecsize + i] = ycoord[i];
            zcoord[enkf_subvecsize + i] = zcoord[i];
         };
      }

      /* store local dimensions for later use */
    nx_local = nx;
    ny_local = ny;
    nz_local = nz;
    origin_local[0] = ix+1;
    origin_local[1] = iy+1;
    origin_local[2] = iz+1;
	}
    //free(xcoord);
    //free(ycoord);
    //free(zcoord);
    //enkf_printvec("info","index", tmpdat);
    //enkf_printvec("info","xcoord", xcoord);
    //enkf_printvec("info","ycoord", ycoord);
    //enkf_printvec("info","zcoord", zcoord);
    //free(tmpdat);
}

void PseudoAdvanceRichards(PFModule *this_module, double start_time, /* Starting time */
double stop_time, /* Stopping time */
PFModule *time_step_control, /* Use this module to control timestep if supplied */
Vector *evap_trans, /* Flux from land surface model */
Vector **pressure_out, /* Output vars */
Vector **porosity_out, Vector **saturation_out);
// gw end

/*-------------------------------------------------------------------------*/
/**
  @author   Wolfgang Kurtz, Guowei He, Mukund Pondkule
  @brief    Initialization of ParFlow for TSMP-PDAF.
  @param    ac    Command line input (for amps).
  @param    *av   Command line input (for amps).
  @param    *input_file   Input file name

  1. First initialization similar to wrf_parflow.c
  2. read time-invariant ET file, if applicable
  3. kuw: create pf vector for printing results to pfb files
  4. create pf vector for printing 2D data
  5. read in mask file (ascii) for overland flow masking
 */
/*--------------------------------------------------------------------------*/
void enkfparflowinit(int ac, char *av[], char *input_file) {

  Grid *grid;

  char *filename = input_file;
  MPI_Comm pfcomm;

  if (screen_wrapper > 1) {
    printf("TSMP-PDAF-WRAPPER mype(w)=%d: enkfparflowinit filename = %s\n", mype_world, filename);
  }

  //L1P_SetStreamPolicy(L1P_stream_optimistic);
  //L1P_SetStreamPolicy(L1P_stream_confirmed);
  //L1P_SetStreamPolicy(L1P_stream_confirmed_or_dcbt);
  //L1P_SetStreamPolicy(L1P_stream_disable);


#ifdef PARFLOW_STAND_ALONE
  pfcomm = MPI_Comm_f2c(comm_model_pdaf);
#endif

  /* BEGINNING: wrf_parflow related part */
  /* ----------------------------------- */

  /*-----------------------------------------------------------------------
   * Initialize AMPS from existing MPI state
   *-----------------------------------------------------------------------*/
#ifdef COUP_OAS_PFL
  if (amps_Init(&ac, &av))
    {
#else
      // Parflow stand alone. No need to guard becasue CLM stand alone should not compile this file.
  if (amps_EmbeddedInitComm(pfcomm))
    {
#endif
      amps_Printf("Error: amps_EmbeddedInit initalization failed\n");
      exit(1);
    }

  /*-----------------------------------------------------------------------
   * Set up globals structure
   *-----------------------------------------------------------------------*/
  NewGlobals(filename);

  /*-----------------------------------------------------------------------
   * Read the Users Input Deck
   *-----------------------------------------------------------------------*/
  amps_ThreadLocal(input_database) = IDB_NewDB(GlobalsInFileName);

  /*-----------------------------------------------------------------------
   * Setup log printing
   *-----------------------------------------------------------------------*/
  NewLogging();

  /*-----------------------------------------------------------------------
   * Setup timing table
   *-----------------------------------------------------------------------*/
  NewTiming();

  /* End of main includes */

  /* Begin of Solver includes */
  GlobalsNumProcsX = GetIntDefault("Process.Topology.P", 1);
  GlobalsNumProcsY = GetIntDefault("Process.Topology.Q", 1);
  GlobalsNumProcsZ = GetIntDefault("Process.Topology.R", 1);

  GlobalsNumProcs = amps_Size(amps_CommWorld);

  GlobalsBackground = ReadBackground();

  GlobalsUserGrid = ReadUserGrid();

  SetBackgroundBounds(GlobalsBackground, GlobalsUserGrid);

  GlobalsMaxRefLevel = 0;

  amps_ThreadLocal(Solver_module) = PFModuleNewModuleType(
							  SolverImpesNewPublicXtraInvoke, SolverRichards, ("Solver"));

  amps_ThreadLocal(solver) = PFModuleNewInstance(
						 amps_ThreadLocal(Solver_module), ());
  /* End of solver includes */

  SetupRichards(amps_ThreadLocal(solver));

  /* Create the flow grid */
  grid = CreateGrid(GlobalsUserGrid);

  /* Create the PF vector holding flux */
  amps_ThreadLocal(evap_trans) = NewVectorType(grid, 1, 1, vector_cell_centered);
  InitVectorAll(amps_ThreadLocal(evap_trans), 0.0);

  /* END: wrf_parflow related part */
  /* ----------------------------- */

  /* read time-invariant ET file, if applicable */
  int etfile = GetEvapTransFile(amps_ThreadLocal(solver));
  if (etfile) {
    char *etfilename = GetEvapTransFilename(amps_ThreadLocal(solver));
    char  filename[256];
    sprintf(filename, "%s", etfilename);
    ReadPFBinary( filename, evap_trans );
    VectorUpdateCommHandle *handle;
    handle = InitVectorUpdate(evap_trans, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
  }

  /* kuw: create pf vector for printing results to pfb files */
  amps_ThreadLocal(vdummy_3d) = NewVectorType(grid, 1, 1, vector_cell_centered);
  InitVectorAll(amps_ThreadLocal(vdummy_3d), 0.0);
  enkf_subvecsize = enkf_getsubvectorsize(grid);
  if(screen_wrapper > 2 && task_id == 1) {
    printf("TSMP-PDAF-WRAPPER mype(w)=%d: enkf_subvecsize=%d\n", mype_world, enkf_subvecsize);
  }

  /* create pf vector for printing 2D data */
  ProblemData *problem_data = GetProblemDataRichards(solver);
  amps_ThreadLocal(vdummy_2d) = NewVectorType(VectorGrid(ProblemDataMannings(problem_data)),1,1,vector_cell_centered);
  InitVectorAll(amps_ThreadLocal(vdummy_2d),0.0);

  /* read in mask file (ascii) for overland flow masking */
  if(pf_olfmasking == 2){
    FILE *friverid=NULL;
    int i;
    friverid = fopen("river.dat","rb");
    fscanf(friverid,"%d",&nriverid);
    riveridx = (int*) calloc(sizeof(int),nriverid);
    riveridy = (int*) calloc(sizeof(int),nriverid);
    for(i=0;i<nriverid;i++){
      fscanf(friverid,"%d %d",&riveridx[i],&riveridy[i]);
      riveridx[i] = riveridx[i]-1;
      riveridy[i] = riveridy[i]-1;
    }
    fclose(friverid);
  }

}

/*-------------------------------------------------------------------------*/
/**
  @author   Wolfgang Kurtz, Guowei He, Mukund Pondkule
  @brief    Initialization of OASIS/ParFlow for TSMP-PDAF.
  @param    current_time    Command line input (for amps).
  @param    dt   Command line input (for amps).

  1. First initialization similar to wrf_parflow.c (wrfparflowadvance_)
  2. Set problem and pressure_in
  3. Allocate and initialize idx_map_subvec2state
  4. Set statevector-size and allocate ParFlow Subvectors
 */
/*--------------------------------------------------------------------------*/
void parflow_oasis_init(double current_time, double dt) {
  double stop_time = current_time + dt;

  Vector *pressure_out;
  Vector *porosity_out;
  Vector *saturation_out;

  VectorUpdateCommHandle *handle;

  /* BEGINNING: wrf_parflow related part */
  /* ----------------------------------- */

  /*
   * Exchange ghost layer data for the newly set fluxes
   */
  handle = InitVectorUpdate(evap_trans, VectorUpdateAll);
  FinalizeVectorUpdate(handle);

  PFModule *time_step_control;

  time_step_control = NewPFModule((void *) SelectTimeStep,
				  (void *) SelectTimeStepInitInstanceXtra,
				  (void *) SelectTimeStepFreeInstanceXtra,
				  (void *) SelectTimeStepNewPublicXtra,
				  (void *) SelectTimeStepFreePublicXtra,
				  (void *) SelectTimeStepSizeOfTempData, NULL, NULL);

  ThisPFModule = time_step_control;
  SelectTimeStepNewPublicXtra();
  ThisPFModule = NULL;

  PFModule *time_step_control_instance = PFModuleNewInstance(
							     time_step_control, ());

  // gw init the OAS, but weird implicit declaration..
  PseudoAdvanceRichards(amps_ThreadLocal(solver), current_time, stop_time,
			time_step_control_instance, amps_ThreadLocal(evap_trans),
			&pressure_out, &porosity_out, &saturation_out);

  PFModuleFreeInstance(time_step_control_instance);
  PFModuleFreeModule(time_step_control);

  /* END: wrf_parflow related part */
  /* ----------------------------- */

  // gw: init idx as well
  Problem *problem = GetProblemRichards(solver);
  Vector *pressure_in = GetPressureRichards(solver);

  /* Allocate and initialize idx_map_subvec2state */
  idx_map_subvec2state   = (int *)   malloc(enkf_subvecsize * sizeof(int));
  init_idx_map_subvec2state(pressure_in);

  /* Set statevector-size and allocate ParFlow Subvectors */
  pf_statevecsize = enkf_subvecsize;
  if(pf_updateflag == 3) pf_statevecsize = pf_statevecsize * 2;
#ifdef FOR2131
  if(pf_updateflag == 2) pf_statevecsize = pf_statevecsize * 2;
#endif
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
  subvec_param_mean      = (double*) calloc(pf_paramvecsize,sizeof(double));
  subvec_param_sd        = (double*) calloc(pf_paramvecsize,sizeof(double));
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
}

/*-------------------------------------------------------------------------*/
/**
  @author   Wolfgang Kurtz, Guowei He, Mukund Pondkule
  @brief    Integration of ParFlow from `current_time` by tim `dt`.
  @param    current_time  Starting time of the simulation.
  @param    dt   Time difference for simulation.

  Integration of ParFlow similar to `wrf_parflow.c` (`wrfparflowadvance_`)

  The `wrf_parflow` related part is essentially the first part of
  `wrfparflowadvance_` without getting `problem_data` and without
  initializing a `PFModule *time_step_control`. Instead there is
  `NULL` given as time step control input (time steps should be
  fixed).

  Afterwards different routines generate a PDAF-statevector from
  simulation outputs. The routines are chosen by flags `pf_updateflag`
  and `pf_paramupdate`.

 */
/*--------------------------------------------------------------------------*/
void enkfparflowadvance(int tcycle, double current_time, double dt)

{
	int i,j;

	/* BEGINNING: wrf_parflow related part */
	/* ----------------------------------- */
	double stop_time = current_time + dt;

	Vector *pressure_out;
	Vector *porosity_out;
	Vector *saturation_out;

	VectorUpdateCommHandle *handle;

	handle = InitVectorUpdate(evap_trans, VectorUpdateAll);
	FinalizeVectorUpdate(handle);

	AdvanceRichards(amps_ThreadLocal(solver), current_time, stop_time, NULL, amps_ThreadLocal(evap_trans), &pressure_out, &porosity_out, &saturation_out);

	handle = InitVectorUpdate(pressure_out, VectorUpdateAll);
	FinalizeVectorUpdate(handle);
	handle = InitVectorUpdate(porosity_out, VectorUpdateAll);
	FinalizeVectorUpdate(handle);
	handle = InitVectorUpdate(saturation_out, VectorUpdateAll);
	FinalizeVectorUpdate(handle);
	/* END: wrf_parflow related part */
	/* ----------------------------- */

	/* create state vector: pressure */
	if(pf_updateflag == 1) {
  	  PF2ENKF(pressure_out, subvec_p);
  	  for(i=0;i<enkf_subvecsize;i++){
             pf_statevec[i] = subvec_p[i];
          }

          /* masking option using saturated cells only */
          if(pf_gwmasking == 1){
            for(i=0;i<enkf_subvecsize;i++){
              subvec_gwind[i] = 1.0;
              if(subvec_p[i]< 0.0) subvec_gwind[i] = 0.0;
            }
          }

          /* masking option using mixed state vector */
          if(pf_gwmasking == 2){
            int no_obs,haveobs,tmpidx;
            MPI_Comm comm_couple_c = MPI_Comm_f2c(comm_couple);
            PF2ENKF(saturation_out, subvec_sat);
  	    PF2ENKF(porosity_out, subvec_porosity);
            MPI_Allreduce(subvec_sat,subvec_mean,enkf_subvecsize,MPI_DOUBLE,MPI_SUM,comm_couple_c);
            for(i=0;i<enkf_subvecsize;i++){
              subvec_gwind[i] = 1.0;
              if(subvec_mean[i]< (double)nreal){
                subvec_gwind[i] = 0.0;
                pf_statevec[i] = subvec_sat[i] * subvec_porosity[i];
              }
            }
            if(task_id == 1 && pf_printgwmask == 1) enkf_printstatistics_pfb(subvec_gwind,"gwind",tstartcycle + stat_dumpoffset,outdir,3);
            get_obsindex_currentobsfile(&no_obs);

            for(i=0;i<no_obs;i++){
              haveobs=0;
              if((xidx_obs[i]>=origin_local[0]) && (xidx_obs[i]<(origin_local[0]+nx_local))){
                if((yidx_obs[i]>=origin_local[1]) && (yidx_obs[i]<(origin_local[1]+ny_local))){
                  if((zidx_obs[i]>=origin_local[2]) && (zidx_obs[i]<(origin_local[2]+nz_local))){
                    haveobs=1;
                  }
                }
              }
              if(haveobs && ind_obs[i]==1){
                for(j=0;j<=(zidx_obs[i]-origin_local[2]);j++){
                  tmpidx = nx_local*ny_local*j + nx_local*(yidx_obs[i]-origin_local[1]) + (xidx_obs[i]-origin_local[0]);
                  subvec_gwind[tmpidx] = 1.0;
                  pf_statevec[tmpidx] = subvec_p[tmpidx];
                }
              }
            }
            if(task_id == 1 && pf_printgwmask == 1) enkf_printstatistics_pfb(subvec_gwind,"gwind_corrected",tstartcycle + stat_dumpoffset,outdir,3);
            clean_obs_pf();
          }
        }

	/* create state vector: swc */
	if(pf_updateflag == 2){
	  PF2ENKF(saturation_out, subvec_sat);
	  PF2ENKF(porosity_out, subvec_porosity);
	  for(i=0;i<enkf_subvecsize;i++) {
            pf_statevec[i] = subvec_sat[i] * subvec_porosity[i];
          }
#ifdef FOR2131
          for(i=enkf_subvecsize,j=0;i<(2*enkf_subvecsize);i++,j++){
               pf_statevec[i] = subvec_p[j];
           }
#endif
	  /* hcp CRNS begins */
	  double dz_glob=GetDouble("ComputationalGrid.DZ");  //hcp if crns update
	  int isc;
	  soilay = (double *) malloc(nz_glob * sizeof(double));
	  char key[IDB_MAX_KEY_LEN];
	  for (isc = 0; isc < nz_glob; isc++) {
	    sprintf(key, "Cell.%d.dzScale.Value", isc);
	    soilay[isc] = GetDouble(key);
	    soilay[isc] *= dz_glob;
	  }
	  /* hcp CRNS ends */
	}

	/* create state vector: joint swc + pressure */
        if(pf_updateflag == 3){
          PF2ENKF(pressure_out, subvec_p);
          PF2ENKF(saturation_out, subvec_sat);
          PF2ENKF(porosity_out, subvec_porosity);
	  for(i=0;i<enkf_subvecsize;i++){
            pf_statevec[i] = subvec_sat[i] * subvec_porosity[i];
          }
          for(i=enkf_subvecsize,j=0;i<(2*enkf_subvecsize);i++,j++){
            pf_statevec[i] = subvec_p[j];
          }
        }

	/* append hydraulic conductivity to state vector */
        if(pf_paramupdate == 1){
           ProblemData *problem_data = GetProblemDataRichards(solver);
           Vector      *perm_xx = ProblemDataPermeabilityX(problem_data);
           handle = InitVectorUpdate(perm_xx, VectorUpdateAll);
           FinalizeVectorUpdate(handle);
           PF2ENKF(perm_xx,subvec_param);
           for(i=(pf_statevecsize-pf_paramvecsize),j=0;i<pf_statevecsize;i++,j++){
             pf_statevec[i] = log10(subvec_param[j]);
           }
        }

	/* append Mannings coefficient to state vector */
        if(pf_paramupdate == 2){
           ProblemData *problem_data = GetProblemDataRichards(solver);
           Vector       *mannings    = ProblemDataMannings(problem_data);
           handle = InitVectorUpdate(mannings, VectorUpdateAll);
           FinalizeVectorUpdate(handle);
           PF2ENKF(mannings,subvec_param);
           for(i=(pf_statevecsize-pf_paramvecsize),j=0;i<pf_statevecsize;i++,j++){
              pf_statevec[i] = log10(subvec_param[j]);
           }
        }

        /* append Porosity to state vector */
        if(pf_paramupdate == 3){
           ProblemData *problem_data = GetProblemDataRichards(solver);
           Vector       *porosity    = ProblemDataPorosity(problem_data);
           handle = InitVectorUpdate(porosity, VectorUpdateAll);
           FinalizeVectorUpdate(handle);
           PF2ENKF(porosity, subvec_param);
           for(i=(pf_statevecsize-pf_paramvecsize),j=0;i<pf_statevecsize;i++,j++){
               pf_statevec[i] = subvec_param[j];
           }
        }

        /* append van Genuchten parameters to state vector */
        if(pf_paramupdate == 4){
           PFModule *relPerm = GetPhaseRelPerm(solver);
           Vector *alpha = PhaseRelPermGetAlpha(relPerm);
           Vector *n = PhaseRelPermGetN(relPerm);
           handle = InitVectorUpdate(alpha, VectorUpdateAll);
           FinalizeVectorUpdate(handle);
           handle = InitVectorUpdate(n, VectorUpdateAll);
           FinalizeVectorUpdate(handle);
           PF2ENKF_2P(alpha, n, subvec_param);
           for(i=(pf_statevecsize-pf_paramvecsize),j=0;i<pf_statevecsize;i++,j++){
               if((j%2)==0){
                   pf_statevec[i] = log(subvec_param[j]);
              }else{
                   pf_statevec[i] = subvec_param[j];
              }
           }
          }

        /* append hydraulic conductivity and porosity to state vector */
        if(pf_paramupdate == 5){
            ProblemData *problem_data = GetProblemDataRichards(solver);
            Vector      *perm_xx = ProblemDataPermeabilityX(problem_data);
            Vector      *porosity    = ProblemDataPorosity(problem_data);
            handle = InitVectorUpdate(perm_xx, VectorUpdateAll);
            FinalizeVectorUpdate(handle);
            handle = InitVectorUpdate(porosity, VectorUpdateAll);
            FinalizeVectorUpdate(handle);
            PF2ENKF_2P(perm_xx,porosity,subvec_param);
            for(i=(pf_statevecsize-pf_paramvecsize),j=0;i<pf_statevecsize;i++,j++){
                    if((j%2)==0){
                        pf_statevec[i] = log10(subvec_param[j]);
                    }else{
                        pf_statevec[i] = subvec_param[j];
                    }
            }
        }

        /* append hydraulic conductivity and van Genuchten parameters to state vector */
        if(pf_paramupdate == 6){
            ProblemData *problem_data = GetProblemDataRichards(solver);
            PFModule    *relPerm = GetPhaseRelPerm(solver);
            Vector      *perm_xx = ProblemDataPermeabilityX(problem_data);
            Vector      *alpha = PhaseRelPermGetAlpha(relPerm);
            Vector      *n     = PhaseRelPermGetN(relPerm);
            handle = InitVectorUpdate(perm_xx, VectorUpdateAll);
            FinalizeVectorUpdate(handle);
            handle = InitVectorUpdate(alpha, VectorUpdateAll);
            FinalizeVectorUpdate(handle);
            handle = InitVectorUpdate(n, VectorUpdateAll);
            FinalizeVectorUpdate(handle);
            PF2ENKF_3P(perm_xx,alpha,n,subvec_param);
            for(i=(pf_statevecsize-pf_paramvecsize),j=0;i<pf_statevecsize;i=i+3,j=j+3){
                    pf_statevec[i] = log10(subvec_param[j]);
                    pf_statevec[i+1] = log(subvec_param[j+1]);
                    pf_statevec[i+2] = subvec_param[j+2];
            }
        }

        /* append porosity and van Genuchten parameters to state vector */
        if(pf_paramupdate == 7){
            ProblemData *problem_data = GetProblemDataRichards(solver);
            PFModule    *relPerm = GetPhaseRelPerm(solver);
            Vector      *porosity    = ProblemDataPorosity(problem_data);
            Vector      *alpha = PhaseRelPermGetAlpha(relPerm);
            Vector      *n     = PhaseRelPermGetN(relPerm);
            handle = InitVectorUpdate(porosity, VectorUpdateAll);
            FinalizeVectorUpdate(handle);
            handle = InitVectorUpdate(alpha, VectorUpdateAll);
            FinalizeVectorUpdate(handle);
            handle = InitVectorUpdate(n, VectorUpdateAll);
            FinalizeVectorUpdate(handle);
            PF2ENKF_3P(porosity,alpha,n,subvec_param);
            for(i=(pf_statevecsize-pf_paramvecsize),j=0;i<pf_statevecsize;i=i+3,j=j+3){
                    pf_statevec[i] = subvec_param[j];
                    pf_statevec[i+1] = log(subvec_param[j+1]);
                    pf_statevec[i+2] = subvec_param[j+2];
            }
        }

        /* append hydraulic conductivity, porosity and van Genuchten parameters to state vector */
        if(pf_paramupdate == 8){
            ProblemData *problem_data = GetProblemDataRichards(solver);
            PFModule    *relPerm = GetPhaseRelPerm(solver);
            Vector      *perm_xx = ProblemDataPermeabilityX(problem_data);
            Vector      *porosity    = ProblemDataPorosity(problem_data);
            Vector      *alpha = PhaseRelPermGetAlpha(relPerm);
            Vector      *n     = PhaseRelPermGetN(relPerm);
            handle = InitVectorUpdate(perm_xx, VectorUpdateAll);
            FinalizeVectorUpdate(handle);
            handle = InitVectorUpdate(porosity, VectorUpdateAll);
            FinalizeVectorUpdate(handle);
            handle = InitVectorUpdate(alpha, VectorUpdateAll);
            FinalizeVectorUpdate(handle);
            handle = InitVectorUpdate(n, VectorUpdateAll);
            FinalizeVectorUpdate(handle);
            PF2ENKF_4P(perm_xx,porosity,alpha,n,subvec_param);
            for(i=(pf_statevecsize-pf_paramvecsize),j=0;i<pf_statevecsize;i=i+4,j=j+4){
                    pf_statevec[i] = log10(subvec_param[j]);
                    pf_statevec[i+1] = subvec_param[j+1];
                    pf_statevec[i+2] = log(subvec_param[j+2]);
                    pf_statevec[i+3] = subvec_param[j+3];
            }
        }

}

void enkfparflowfinalize() {

	// gw: free assimilation data structures
	free(idx_map_subvec2state);
	free(xcoord);
	free(ycoord);
	free(zcoord);
	/* hcp CRNS begins */
	if(pf_updateflag == 2){
	  free(soilay);
	}
	/* hcp CRNS ends */

	fflush(NULL);
	LogGlobals();
	PrintTiming();
	FreeLogging();
	FreeTiming();
	FreeGlobals();
	amps_Finalize();
}

int enkf_getsubvectorsize(Grid *grid) {
	int sg;
	int out = 0;
	ForSubgridI(sg, GridSubgrids(grid))
	{
		Subgrid *subgrid = GridSubgrid(grid, sg);
		int nx = SubgridNX(subgrid);
		int ny = SubgridNY(subgrid);
		int nz = SubgridNZ(subgrid);
		out += nx * ny * nz; // kuw: correct?
	}
	return (out);
}

void PF2ENKF(Vector *pf_vector, double *enkf_subvec) {

	Grid *grid = VectorGrid(pf_vector);
	int sg;

	ForSubgridI(sg, GridSubgrids(grid))
	{
		Subgrid *subgrid = GridSubgrid(grid, sg);

		int ix = SubgridIX(subgrid);
		int iy = SubgridIY(subgrid);
		int iz = SubgridIZ(subgrid);

		int nx = SubgridNX(subgrid);
		int ny = SubgridNY(subgrid);
		int nz = SubgridNZ(subgrid);


		Subvector *subvector = VectorSubvector(pf_vector, sg);
		double *subvector_data = SubvectorData(subvector);

		int i, j, k;
		int counter = 0;

		for (k = iz; k < iz + nz; k++) {
			for (j = iy; j < iy + ny; j++) {
				for (i = ix; i < ix + nx; i++) {
					int pf_index = SubvectorEltIndex(subvector, i, j, k);
					enkf_subvec[counter] = (subvector_data[pf_index]);
					counter++;
				}
			}
		}
	}
}

void PF2ENKF_2P(Vector *p1_vector, Vector *p2_vector, double *enkf_subvec) {

	Grid *grid = VectorGrid(p1_vector);
	int sg;

	ForSubgridI(sg, GridSubgrids(grid))
	{
		Subgrid *subgrid = GridSubgrid(grid, sg);

		int ix = SubgridIX(subgrid);
		int iy = SubgridIY(subgrid);
		int iz = SubgridIZ(subgrid);

		int nx = SubgridNX(subgrid);
		int ny = SubgridNY(subgrid);
		int nz = SubgridNZ(subgrid);


		Subvector *subvector = VectorSubvector(p1_vector, sg);
		double *subvector_data = SubvectorData(subvector);

                Subvector *subvector_2 = VectorSubvector(p2_vector, sg);
		double *subvector_2_data = SubvectorData(subvector_2);

		int i, j, k;
		int counter = 0;

		for (k = iz; k < iz + nz; k++) {
			for (j = iy; j < iy + ny; j++) {
				for (i = ix; i < ix + nx; i++) {
					int pf_index = SubvectorEltIndex(subvector, i, j, k);
                                        int pf_index_2 = SubvectorEltIndex(subvector_2, i, j, k);
					enkf_subvec[counter] = (subvector_data[pf_index]);
					counter++;
                                        enkf_subvec[counter] = (subvector_2_data[pf_index_2]);
					counter++;
				}
			}
		}
	}

}

void PF2ENKF_3P(Vector *p1_vector, Vector *p2_vector, Vector *p3_vector, double *enkf_subvec) {

	Grid *grid = VectorGrid(p1_vector);
	int sg;

	ForSubgridI(sg, GridSubgrids(grid))
	{
		Subgrid *subgrid = GridSubgrid(grid, sg);

		int ix = SubgridIX(subgrid);
		int iy = SubgridIY(subgrid);
		int iz = SubgridIZ(subgrid);

		int nx = SubgridNX(subgrid);
		int ny = SubgridNY(subgrid);
		int nz = SubgridNZ(subgrid);


		Subvector *subvector = VectorSubvector(p1_vector, sg);
		double *subvector_data = SubvectorData(subvector);

                Subvector *subvector_2 = VectorSubvector(p2_vector, sg);
		double *subvector_2_data = SubvectorData(subvector_2);

                Subvector *subvector_3 = VectorSubvector(p3_vector, sg);
		double *subvector_3_data = SubvectorData(subvector_3);

		int i, j, k;
		int counter = 0;

		for (k = iz; k < iz + nz; k++) {
			for (j = iy; j < iy + ny; j++) {
				for (i = ix; i < ix + nx; i++) {
					int pf_index = SubvectorEltIndex(subvector, i, j, k);
                                        int pf_index_2 = SubvectorEltIndex(subvector_2, i, j, k);
                                        int pf_index_3 = SubvectorEltIndex(subvector_3, i, j, k);
					enkf_subvec[counter] = (subvector_data[pf_index]);
					counter++;
                                        enkf_subvec[counter] = (subvector_2_data[pf_index_2]);
					counter++;
                                        enkf_subvec[counter] = (subvector_3_data[pf_index_3]);
					counter++;
				}
			}
		}
	}

}

void PF2ENKF_4P(Vector *p1_vector, Vector *p2_vector, Vector *p3_vector, Vector *p4_vector, double *enkf_subvec) {

	Grid *grid = VectorGrid(p1_vector);
	int sg;

	ForSubgridI(sg, GridSubgrids(grid))
	{
		Subgrid *subgrid = GridSubgrid(grid, sg);

		int ix = SubgridIX(subgrid);
		int iy = SubgridIY(subgrid);
		int iz = SubgridIZ(subgrid);

		int nx = SubgridNX(subgrid);
		int ny = SubgridNY(subgrid);
		int nz = SubgridNZ(subgrid);


		Subvector *subvector = VectorSubvector(p1_vector, sg);
		double *subvector_data = SubvectorData(subvector);

                Subvector *subvector_2 = VectorSubvector(p2_vector, sg);
		double *subvector_2_data = SubvectorData(subvector_2);

                Subvector *subvector_3 = VectorSubvector(p3_vector, sg);
		double *subvector_3_data = SubvectorData(subvector_3);

                Subvector *subvector_4 = VectorSubvector(p4_vector, sg);
		double *subvector_4_data = SubvectorData(subvector_4);

		int i, j, k;
		int counter = 0;

		for (k = iz; k < iz + nz; k++) {
			for (j = iy; j < iy + ny; j++) {
				for (i = ix; i < ix + nx; i++) {
					int pf_index = SubvectorEltIndex(subvector, i, j, k);
                                        int pf_index_2 = SubvectorEltIndex(subvector_2, i, j, k);
                                        int pf_index_3 = SubvectorEltIndex(subvector_3, i, j, k);
                                        int pf_index_4 = SubvectorEltIndex(subvector_4, i, j, k);
					enkf_subvec[counter] = (subvector_data[pf_index]);
					counter++;
                                        enkf_subvec[counter] = (subvector_2_data[pf_index_2]);
					counter++;
                                        enkf_subvec[counter] = (subvector_3_data[pf_index_3]);
					counter++;
                                        enkf_subvec[counter] = (subvector_4_data[pf_index_4]);
					counter++;
				}
			}
		}
	}

}

void ENKF2PF(Vector *pf_vector, double *enkf_subvec) {

	Grid *grid = VectorGrid(amps_ThreadLocal(pf_vector));
	int sg;

	ForSubgridI(sg, GridSubgrids(grid))
	{
		Subgrid *subgrid = GridSubgrid(grid, sg);

		int ix = SubgridIX(subgrid);
		int iy = SubgridIY(subgrid);
		int iz = SubgridIZ(subgrid);

		int nx = SubgridNX(subgrid);
		int ny = SubgridNY(subgrid);
		int nz = SubgridNZ(subgrid);

		Subvector *subvector = VectorSubvector(pf_vector, sg);
		double *subvector_data = SubvectorData(subvector);

		int i, j, k;
		int counter = 0;

		for (k = iz; k < iz + nz; k++) {
			for (j = iy; j < iy + ny; j++) {
				for (i = ix; i < ix + nx; i++) {
					int pf_index = SubvectorEltIndex(subvector, i, j, k);
					(subvector_data[pf_index]) = enkf_subvec[counter];
					counter++;
				}
			}
		}
	}
}

void ENKF2PF_2P(Vector *p1_vector, Vector *p2_vector, double *enkf_subvec) {

	Grid *grid = VectorGrid(amps_ThreadLocal(p1_vector));
	int sg;

	ForSubgridI(sg, GridSubgrids(grid))
	{
		Subgrid *subgrid = GridSubgrid(grid, sg);

		int ix = SubgridIX(subgrid);
		int iy = SubgridIY(subgrid);
		int iz = SubgridIZ(subgrid);

		int nx = SubgridNX(subgrid);
		int ny = SubgridNY(subgrid);
		int nz = SubgridNZ(subgrid);

		Subvector *subvector = VectorSubvector(p1_vector, sg);
		double *subvector_data = SubvectorData(subvector);

                Subvector *subvector_2 = VectorSubvector(p2_vector, sg);
		double *subvector_2_data = SubvectorData(subvector_2);

		int i, j, k;
		int counter = 0;

		for (k = iz; k < iz + nz; k++) {
			for (j = iy; j < iy + ny; j++) {
				for (i = ix; i < ix + nx; i++) {
					int pf_index = SubvectorEltIndex(subvector, i, j, k);
                                        int pf_index_2 = SubvectorEltIndex(subvector_2, i, j, k);
					(subvector_data[pf_index]) = enkf_subvec[counter];
					counter++;
                                        (subvector_2_data[pf_index_2]) = enkf_subvec[counter];
					counter++;
				}
			}
		}
	}
}

void ENKF2PF_3P(Vector *p1_vector, Vector *p2_vector, Vector *p3_vector, double *enkf_subvec) {

	Grid *grid = VectorGrid(amps_ThreadLocal(p1_vector));
	int sg;

	ForSubgridI(sg, GridSubgrids(grid))
	{
		Subgrid *subgrid = GridSubgrid(grid, sg);

		int ix = SubgridIX(subgrid);
		int iy = SubgridIY(subgrid);
		int iz = SubgridIZ(subgrid);

		int nx = SubgridNX(subgrid);
		int ny = SubgridNY(subgrid);
		int nz = SubgridNZ(subgrid);

		Subvector *subvector = VectorSubvector(p1_vector, sg);
		double *subvector_data = SubvectorData(subvector);

                Subvector *subvector_2 = VectorSubvector(p2_vector, sg);
		double *subvector_2_data = SubvectorData(subvector_2);

                Subvector *subvector_3 = VectorSubvector(p3_vector, sg);
		double *subvector_3_data = SubvectorData(subvector_3);

		int i, j, k;
		int counter = 0;

		for (k = iz; k < iz + nz; k++) {
			for (j = iy; j < iy + ny; j++) {
				for (i = ix; i < ix + nx; i++) {
					int pf_index = SubvectorEltIndex(subvector, i, j, k);
                                        int pf_index_2 = SubvectorEltIndex(subvector_2, i, j, k);
                                        int pf_index_3 = SubvectorEltIndex(subvector_3, i, j, k);
					(subvector_data[pf_index]) = enkf_subvec[counter];
					counter++;
                                        (subvector_2_data[pf_index_2]) = enkf_subvec[counter];
					counter++;
                                        (subvector_3_data[pf_index_3]) = enkf_subvec[counter];
					counter++;
				}
			}
		}
	}
}

void ENKF2PF_4P(Vector *p1_vector, Vector *p2_vector, Vector *p3_vector, Vector *p4_vector, double *enkf_subvec) {

	Grid *grid = VectorGrid(amps_ThreadLocal(p1_vector));
	int sg;

	ForSubgridI(sg, GridSubgrids(grid))
	{
		Subgrid *subgrid = GridSubgrid(grid, sg);

		int ix = SubgridIX(subgrid);
		int iy = SubgridIY(subgrid);
		int iz = SubgridIZ(subgrid);

		int nx = SubgridNX(subgrid);
		int ny = SubgridNY(subgrid);
		int nz = SubgridNZ(subgrid);

		Subvector *subvector = VectorSubvector(p1_vector, sg);
		double *subvector_data = SubvectorData(subvector);

                Subvector *subvector_2 = VectorSubvector(p2_vector, sg);
		double *subvector_2_data = SubvectorData(subvector_2);

                Subvector *subvector_3 = VectorSubvector(p3_vector, sg);
		double *subvector_3_data = SubvectorData(subvector_3);

                Subvector *subvector_4 = VectorSubvector(p4_vector, sg);
		double *subvector_4_data = SubvectorData(subvector_4);

		int i, j, k;
		int counter = 0;

		for (k = iz; k < iz + nz; k++) {
			for (j = iy; j < iy + ny; j++) {
				for (i = ix; i < ix + nx; i++) {
					int pf_index = SubvectorEltIndex(subvector, i, j, k);
                                        int pf_index_2 = SubvectorEltIndex(subvector_2, i, j, k);
                                        int pf_index_3 = SubvectorEltIndex(subvector_3, i, j, k);
                                        int pf_index_4 = SubvectorEltIndex(subvector_4, i, j, k);
					(subvector_data[pf_index]) = enkf_subvec[counter];
					counter++;
                                        (subvector_2_data[pf_index_2]) = enkf_subvec[counter];
					counter++;
                                        (subvector_3_data[pf_index_3]) = enkf_subvec[counter];
					counter++;
                                        (subvector_4_data[pf_index_4]) = enkf_subvec[counter];
					counter++;
				}
			}
		}
	}
}

void ENKF2PF_masked(Vector *pf_vector, double *enkf_subvec, double *mask) {

	Grid *grid = VectorGrid(pf_vector);
	int sg;

	ForSubgridI(sg, GridSubgrids(grid))
	{
		Subgrid *subgrid = GridSubgrid(grid, sg);

		int ix = SubgridIX(subgrid);
		int iy = SubgridIY(subgrid);
		int iz = SubgridIZ(subgrid);

		int nx = SubgridNX(subgrid);
		int ny = SubgridNY(subgrid);
		int nz = SubgridNZ(subgrid);

		Subvector *subvector = VectorSubvector(pf_vector, sg);
		double *subvector_data = SubvectorData(subvector);

		int i, j, k;
		int counter = 0;

		for (k = iz; k < iz + nz; k++) {
			for (j = iy; j < iy + ny; j++) {
				for (i = ix; i < ix + nx; i++) {
					int pf_index = SubvectorEltIndex(subvector, i, j, k);
					subvector_data[pf_index] = mask[counter]*enkf_subvec[counter] + (1.0-mask[counter])*subvector_data[pf_index];
					counter++;
				}
			}
		}
	}
}

void enkf_printvec(char *pre, char *suff, double *data, int dim) {
  Vector *v=NULL;
  if(dim==2){
    v = vdummy_2d;
  }else{
    v = vdummy_3d;
  }
  Grid *grid = VectorGrid(v);
  int sg;

  ForSubgridI(sg, GridSubgrids(grid))
  {
    Subgrid *subgrid = GridSubgrid(grid, sg);
    int ix = SubgridIX(subgrid);
    int iy = SubgridIY(subgrid);
    int iz = SubgridIZ(subgrid);

    int nx = SubgridNX(subgrid);
    int ny = SubgridNY(subgrid);
    int nz = SubgridNZ(subgrid);

    Subvector *subvector = VectorSubvector(v, sg);

    double *subvector_data = SubvectorData(subvector);
    int     i, j, k;
    int     counter = 0;

    for (k = iz; k < iz + nz; k++) {
      for (j = iy; j < iy + ny; j++) {
    	for (i = ix; i < ix + nx; i++) {
    	  int pf_index = SubvectorEltIndex(subvector, i, j, k);
    	  subvector_data[pf_index] = data[counter];
    	  counter++;
    	}
      }
    }

  }

  WritePFBinary(pre, suff, v);
}

void enkf_printmannings(char *pre, char *suff){
    ProblemData *problem_data = GetProblemDataRichards(solver);
    WritePFBinary(pre,suff, ProblemDataMannings(problem_data));
}


void update_parflow (int do_pupd) {
  int i,j;
  VectorUpdateCommHandle *handle;

  if(pf_paramupdate == 3 && do_pupd){
    ProblemData *problem_data = GetProblemDataRichards(solver);
    Vector      *porosity     = ProblemDataPorosity(problem_data);
    int nshift = 0;
    if(pf_updateflag == 3 || pf_updateflag == 2){
      nshift = 2*enkf_subvecsize;
    }else{
      nshift = enkf_subvecsize;
    }
    /* update porosity */
    for(i=nshift,j=0;i<(nshift+pf_paramvecsize);i++,j++){
      subvec_param[j] = pf_statevec[i];
      subvec_porosity[j] = pf_statevec[i];
    }

    ENKF2PF(porosity,subvec_param);
    handle = InitVectorUpdate(porosity,VectorUpdateAll);
    FinalizeVectorUpdate(handle);
  }

  /* hydraulic conductivity and porosity */
  if(pf_paramupdate == 5 && do_pupd){
    ProblemData * problem_data = GetProblemDataRichards(solver);
    Vector            *perm_xx = ProblemDataPermeabilityX(problem_data);
    Vector            *perm_yy = ProblemDataPermeabilityY(problem_data);
    Vector            *perm_zz = ProblemDataPermeabilityZ(problem_data);
    Vector       *porosity    = ProblemDataPorosity(problem_data);
    int nshift = 0;
    if(pf_updateflag == 3 || pf_updateflag == 2){
      nshift = 2*enkf_subvecsize;
    }else{
      nshift = enkf_subvecsize;
    }

    int por_counter = 0;
    for(i=nshift,j=0;i<(nshift+pf_paramvecsize);i++,j++){
        subvec_param[j] = pf_statevec[i];
        if((j%2)!=0){
            subvec_porosity[por_counter] = pf_statevec[i];
            por_counter++;
        }
    }

    ENKF2PF_2P(perm_xx,porosity,subvec_param);
    handle = InitVectorUpdate(perm_xx, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
    handle = InitVectorUpdate(porosity, VectorUpdateAll);
    FinalizeVectorUpdate(handle);

    /* update perm_yy and perm_zz */
    for(i=nshift,j=0;i<(nshift+pf_paramvecsize);i++,j++){
        if((j%2)==0){
            subvec_param[j] = pf_statevec[i] * pf_aniso_perm_y;
            subvec_param[j+1] = pf_statevec[i] * pf_aniso_perm_z;
        }
    }

    ENKF2PF_2P(perm_yy,perm_zz,subvec_param);
    handle = InitVectorUpdate(perm_yy, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
    handle = InitVectorUpdate(perm_zz, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
  }

  /* porosity and van Genuchten parameter */
  if(pf_paramupdate == 7 && do_pupd){
    ProblemData * problem_data = GetProblemDataRichards(solver);
    Vector            *porosity = ProblemDataPorosity(problem_data);
    PFModule *relPerm = GetPhaseRelPerm(solver);
    Vector            *alpha = PhaseRelPermGetAlpha(relPerm);
    Vector            *n = PhaseRelPermGetN(relPerm);
    PFModule *sat     = GetSaturation(solver);
    Vector            *alpha_sat = SaturationGetAlpha(sat);
    Vector            *n_sat = SaturationGetN(sat);
    int nshift = 0;
    if(pf_updateflag == 3 || pf_updateflag == 2){
      nshift = 2*enkf_subvecsize;
    }else{
      nshift = enkf_subvecsize;
    }

    int por_counter = 0;
    for(i=nshift,j=0;i<(nshift+pf_paramvecsize);i++,j++){
        subvec_param[j] = pf_statevec[i];
        if((j%3)==0){
            subvec_porosity[por_counter] = pf_statevec[i];
            por_counter++;
        }
    }

    ENKF2PF_3P(porosity,alpha,n,subvec_param);
    handle = InitVectorUpdate(porosity, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
    handle = InitVectorUpdate(alpha, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
    handle = InitVectorUpdate(n, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
    ENKF2PF_3P(porosity,alpha_sat,n_sat,subvec_param);
    handle = InitVectorUpdate(porosity, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
    handle = InitVectorUpdate(alpha_sat, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
    handle = InitVectorUpdate(n_sat, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
  }

  /* hydraulic conductivity, porosity and van Genuchten parameter */
  if(pf_paramupdate == 8 && do_pupd){
    ProblemData * problem_data = GetProblemDataRichards(solver);
    Vector            *perm_xx = ProblemDataPermeabilityX(problem_data);
    Vector            *perm_yy = ProblemDataPermeabilityY(problem_data);
    Vector            *perm_zz = ProblemDataPermeabilityZ(problem_data);
    Vector            *porosity = ProblemDataPorosity(problem_data);
    PFModule *relPerm = GetPhaseRelPerm(solver);
    Vector            *alpha = PhaseRelPermGetAlpha(relPerm);
    Vector            *n = PhaseRelPermGetN(relPerm);
    PFModule *sat     = GetSaturation(solver);
    Vector            *alpha_sat = SaturationGetAlpha(sat);
    Vector            *n_sat = SaturationGetN(sat);
    int nshift = 0;
    if(pf_updateflag == 3 || pf_updateflag == 2){
      nshift = 2*enkf_subvecsize;
    }else{
      nshift = enkf_subvecsize;
    }

    int por_counter = 0;
    for(i=nshift,j=0;i<(nshift+pf_paramvecsize);i++,j++){
        subvec_param[j] = pf_statevec[i];
        if((j%4)==1){
            subvec_porosity[por_counter] = pf_statevec[i];
            por_counter++;
        }
    }

    ENKF2PF_4P(perm_xx,porosity,alpha,n,subvec_param);
    handle = InitVectorUpdate(perm_xx, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
    handle = InitVectorUpdate(porosity, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
    handle = InitVectorUpdate(alpha, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
    handle = InitVectorUpdate(n, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
    ENKF2PF_4P(perm_xx,porosity,alpha_sat,n_sat,subvec_param);
    handle = InitVectorUpdate(perm_xx, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
    handle = InitVectorUpdate(porosity, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
    handle = InitVectorUpdate(alpha_sat, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
    handle = InitVectorUpdate(n_sat, VectorUpdateAll);
    FinalizeVectorUpdate(handle);

    /* update perm_yy and perm_zz*/
    for(i=nshift,j=0;i<(nshift+pf_paramvecsize);i=i+4,j=j+4){
        subvec_param[j] = pf_statevec[i];
        subvec_param[j+1] = pf_statevec[i] * pf_aniso_perm_y;
        subvec_param[j+2] = pf_statevec[i] * pf_aniso_perm_z;
        subvec_param[j+3] = pf_statevec[i+1];
    }

    ENKF2PF_4P(perm_xx,perm_yy,perm_zz,porosity,subvec_param);
    handle = InitVectorUpdate(perm_xx, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
    handle = InitVectorUpdate(porosity, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
    handle = InitVectorUpdate(perm_yy, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
    handle = InitVectorUpdate(perm_zz, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
  }

  if(pf_olfmasking == 1) mask_overlandcells();
  if(pf_olfmasking == 2) mask_overlandcells_river();

  if(pf_updateflag == 1) {
    Vector *pressure_in = GetPressureRichards(solver);
    //Vector *pressure_in = NULL;
    //pressure_in = GetPressureRichards(amps_ThreadLocal(solver));
    //GetPressureRichards3(amps_ThreadLocal(solver),&pressure_in);
    //pressure_in = (void*) GetPressureRichards4(amps_ThreadLocal(solver));

    /* no groundwater masking */
    if(pf_gwmasking == 0){
      ENKF2PF(pressure_in, pf_statevec);
    }

    /* groundwater masking using saturated cells only */
    if(pf_gwmasking == 1){
      ENKF2PF_masked(pressure_in, pf_statevec,subvec_gwind);
    }

    /* groundwater masking using mixed state vector */
    if(pf_gwmasking == 2){
      Problem      *problem = GetProblemRichards(solver);
      ProblemData  *problem_data = GetProblemDataRichards(solver);
      PFModule     *problem_saturation = ProblemSaturation(problem);
      double       gravity  = ProblemGravity(problem);
      Vector *saturation_in = GetSaturationRichards(solver);
      Vector *density       = GetDensityRichards(solver);
      int saturation_to_pressure_type = 1;

      /* first update swc cells from mixed state vector pf_statevec */
      for(i=0;i<enkf_subvecsize;i++){
        subvec_sat[i] = subvec_gwind[i]*subvec_sat[i] + (1.0-subvec_gwind[i])*pf_statevec[i]/subvec_porosity[i];
      }
      ENKF2PF(saturation_in,subvec_sat);
      global_ptr_this_pf_module = problem_saturation;
      SaturationToPressure(saturation_in,	pressure_in, density, gravity,problem_data, CALCFCN, saturation_to_pressure_type);
      global_ptr_this_pf_module = solver;

      /* second update remaining pressures cells from mixed state vector pf_statevec */
      ENKF2PF_masked(pressure_in,pf_statevec,subvec_gwind);
    }

    /* update ghost cells for pressure */
    handle = InitVectorUpdate(pressure_in, VectorUpdateAll);
    FinalizeVectorUpdate(handle);

  }

  if(pf_updateflag == 2){
    /* write state vector to saturation in parflow */
    Vector * saturation_in = GetSaturationRichards(solver);
    for(i=0;i<enkf_subvecsize;i++){
      pf_statevec[i] = pf_statevec[i] / subvec_porosity[i];
#ifdef FOR2131
      if(pf_statevec[i] > 1.0){
         pf_statevec[i] = 1.0;
      }
#endif
    }
    int saturation_to_pressure_type = 1;
    ENKF2PF(saturation_in, pf_statevec);
    Problem * problem = GetProblemRichards(solver);
    double gravity = ProblemGravity(problem);
    Vector * pressure_in = GetPressureRichards(solver);
    Vector * density = GetDensityRichards(solver);
    ProblemData * problem_data = GetProblemDataRichards(solver);
    PFModule * problem_saturation = ProblemSaturation(problem);
    // convert saturation to pressure
    global_ptr_this_pf_module = problem_saturation;
    SaturationToPressure(saturation_in,	pressure_in, density, gravity, problem_data, CALCFCN, saturation_to_pressure_type);
    global_ptr_this_pf_module = solver;

    PF2ENKF(pressure_in,subvec_p);
    handle = InitVectorUpdate(pressure_in, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
  }

  if(pf_updateflag == 3){
    Vector *pressure_in = GetPressureRichards(solver);

    ENKF2PF(pressure_in,&pf_statevec[enkf_subvecsize]);

    handle = InitVectorUpdate(pressure_in, VectorUpdateAll);
    FinalizeVectorUpdate(handle);

  }

  if(pf_paramupdate == 1 && do_pupd){
    ProblemData * problem_data = GetProblemDataRichards(solver);
    Vector            *perm_xx = ProblemDataPermeabilityX(problem_data);
    Vector            *perm_yy = ProblemDataPermeabilityY(problem_data);
    Vector            *perm_zz = ProblemDataPermeabilityZ(problem_data);
    int nshift = 0;
    if(pf_updateflag == 3){
      nshift = 2*enkf_subvecsize;
    }else{
#ifdef FOR2131
      if(pf_updateflag == 2){
	nshift = 2*enkf_subvecsize;
      }else{
#endif
	nshift = enkf_subvecsize;
#ifdef FOR2131
      }
#endif
    }

    /* update perm_xx */
    for(i=nshift,j=0;i<(nshift+enkf_subvecsize);i++,j++)
      subvec_param[j] = pf_statevec[i];

    if(pf_gwmasking == 0){
      ENKF2PF(perm_xx,subvec_param);
    }
    // hcp gmasking with param
    if(pf_gwmasking == 1){
//      printf("Kxx masked");
      ENKF2PF_masked(perm_xx, subvec_param,subvec_gwind);
    }
    if(pf_gwmasking == 2){
//      printf("Kxx masked");
      ENKF2PF_masked(perm_xx, subvec_param,subvec_gwind);
    }
    // hcp fin
    handle = InitVectorUpdate(perm_xx, VectorUpdateAll);
    FinalizeVectorUpdate(handle);

    /* update perm_yy */
    for(i=nshift,j=0;i<(nshift+enkf_subvecsize);i++,j++)
      subvec_param[j] = pf_statevec[i] * pf_aniso_perm_y;

    if(pf_gwmasking == 0){
      ENKF2PF(perm_yy,subvec_param);
    }
    // hcp gmasking with param
    if(pf_gwmasking == 1){
//      printf("Kyy masked");
      ENKF2PF_masked(perm_yy, subvec_param,subvec_gwind);
    }
    if(pf_gwmasking == 2){
//      printf("Kyy masked");
      ENKF2PF_masked(perm_yy, subvec_param,subvec_gwind);
    }
    // hcp fin
    handle = InitVectorUpdate(perm_yy, VectorUpdateAll);
    FinalizeVectorUpdate(handle);

    /* update perm_zz */
    for(i=nshift,j=0;i<(nshift+enkf_subvecsize);i++,j++)
      subvec_param[j] = pf_statevec[i] * pf_aniso_perm_z;

    if(pf_gwmasking == 0){
      ENKF2PF(perm_zz,subvec_param);
    }
    // hcp gmasking with param
    if(pf_gwmasking == 1){
//      printf("Kzz masked");
      ENKF2PF_masked(perm_zz, subvec_param,subvec_gwind);
    }
    if(pf_gwmasking == 2){
//      printf("Kzz masked");
      ENKF2PF_masked(perm_zz, subvec_param,subvec_gwind);
    }
    // hcp fin
    handle = InitVectorUpdate(perm_zz, VectorUpdateAll);
    FinalizeVectorUpdate(handle);

  }

  if(pf_paramupdate == 2 && do_pupd){
    ProblemData *problem_data = GetProblemDataRichards(solver);
    Vector       *mannings    = ProblemDataMannings(problem_data);
    int nshift = 0;
    if(pf_updateflag == 3){
      nshift = 2*enkf_subvecsize;
    }else{
#ifdef FOR2131
      if(pf_updateflag == 2){
	nshift = 2*enkf_subvecsize;
      }else{
#endif
	nshift = enkf_subvecsize;
#ifdef FOR2131
      }
#endif
    }
    /* update mannings */
    for(i=nshift,j=0;i<(nshift+pf_paramvecsize);i++,j++)
      subvec_param[j] = pf_statevec[i];

    ENKF2PF(mannings,subvec_param);
    handle = InitVectorUpdate(mannings, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
  }

  /* van Genuchten */
  if(pf_paramupdate == 4 && do_pupd){
    PFModule *relPerm = GetPhaseRelPerm(solver);
    Vector            *alpha = PhaseRelPermGetAlpha(relPerm);
    Vector            *n = PhaseRelPermGetN(relPerm);
    PFModule *sat     = GetSaturation(solver);
    Vector            *alpha_sat = SaturationGetAlpha(sat);
    Vector            *n_sat = SaturationGetN(sat);
    int nshift = 0;
    if(pf_updateflag == 3 || pf_updateflag == 2){
      nshift = 2*enkf_subvecsize;
    }else{
      nshift = enkf_subvecsize;
    }

    for(i=nshift,j=0;i<(nshift+pf_paramvecsize);i++,j++)
      subvec_param[j] = pf_statevec[i];

    ENKF2PF_2P(alpha,n,subvec_param);
    handle = InitVectorUpdate(alpha, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
    handle = InitVectorUpdate(n, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
    ENKF2PF_2P(alpha_sat,n_sat,subvec_param);
    handle = InitVectorUpdate(alpha_sat, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
    handle = InitVectorUpdate(n_sat, VectorUpdateAll);
    FinalizeVectorUpdate(handle);

  }

  /* hydraulic conductivity and van Genuchten parameter */
  if(pf_paramupdate == 6 && do_pupd){
    ProblemData * problem_data = GetProblemDataRichards(solver);
    Vector            *perm_xx = ProblemDataPermeabilityX(problem_data);
    Vector            *perm_yy = ProblemDataPermeabilityY(problem_data);
    Vector            *perm_zz = ProblemDataPermeabilityZ(problem_data);
    PFModule *relPerm = GetPhaseRelPerm(solver);
    Vector            *alpha = PhaseRelPermGetAlpha(relPerm);
    Vector            *n = PhaseRelPermGetN(relPerm);
    PFModule *sat     = GetSaturation(solver);
    Vector            *alpha_sat = SaturationGetAlpha(sat);
    Vector            *n_sat = SaturationGetN(sat);
    int nshift = 0;
    if(pf_updateflag == 3 || pf_updateflag == 2){
      nshift = 2*enkf_subvecsize;
    }else{
      nshift = enkf_subvecsize;
    }

    for(i=nshift,j=0;i<(nshift+pf_paramvecsize);i++,j++){
        subvec_param[j] = pf_statevec[i];
    }

    ENKF2PF_3P(perm_xx,alpha,n,subvec_param);
    handle = InitVectorUpdate(perm_xx, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
    handle = InitVectorUpdate(alpha, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
    handle = InitVectorUpdate(n, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
    ENKF2PF_3P(perm_xx,alpha_sat,n_sat,subvec_param);
    handle = InitVectorUpdate(perm_xx, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
    handle = InitVectorUpdate(alpha_sat, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
    handle = InitVectorUpdate(n_sat, VectorUpdateAll);
    FinalizeVectorUpdate(handle);

    /* update perm_yy and perm_zz*/
    for(i=nshift,j=0;i<(nshift+pf_paramvecsize);i=i+3,j=j+3){
        if((j%2)==0){
            subvec_param[j] = pf_statevec[i];
            subvec_param[j+1] = pf_statevec[i] * pf_aniso_perm_y;
            subvec_param[j+2] = pf_statevec[i] * pf_aniso_perm_z;
        }
    }

    ENKF2PF_3P(perm_xx,perm_yy,perm_zz,subvec_param);
    handle = InitVectorUpdate(perm_xx, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
    handle = InitVectorUpdate(perm_yy, VectorUpdateAll);
    FinalizeVectorUpdate(handle);
    handle = InitVectorUpdate(perm_zz, VectorUpdateAll);
    FinalizeVectorUpdate(handle);

  }

}


void mask_overlandcells()
{
  int i,j,k;
  int counter = 0;

  /* fast-forward counter to uppermost model layer */
  counter = nx_local*ny_local*(nz_local-1);

  /* mask updated values in uppermost model layer */
  if(pf_updateflag == 1){
    for(i=0;i<ny_local;i++){
      for(j=0;j<nx_local;j++){
        //if(subvec_p[counter]>0.0) pf_statevec[counter] = subvec_p[counter];
        pf_statevec[counter] = subvec_p[counter];
        counter++;
      }
    }
    if(pf_gwmasking == 2){   //There are overland cells being unsat (by hcp)
      counter = nx_local*ny_local*(nz_local-1);
      for(i=0;i<ny_local;i++){
        for(j=0;j<nx_local;j++){
          //if(subvec_p[counter]>0.0) pf_statevec[counter] = subvec_p[counter];
          if(subvec_gwind[counter] < 0.5){
             pf_statevec[counter] = subvec_sat[counter]*subvec_porosity[counter];
          }
          counter++;
        }
      }
    }
  }
  if(pf_updateflag == 2){
    for(i=0;i<ny_local;i++){
      for(j=0;j<nx_local;j++){
        //if(subvec_p[counter]>0.0) pf_statevec[counter] = subvec_sat[counter]*subvec_porosity[counter];
        pf_statevec[counter] = subvec_sat[counter]*subvec_porosity[counter];
        counter++;
      }
    }
  }
  if(pf_updateflag == 3){
    for(i=0;i<ny_local;i++){
      for(j=0;j<nx_local;j++){
        //if(subvec_p[counter]>0.0) pf_statevec[counter+enkf_subvecsize] = subvec_p[counter];
        pf_statevec[counter+enkf_subvecsize] = subvec_p[counter];
        counter++;
      }
    }
  }

}


void mask_overlandcells_river()
{
  int i,j,idx;

  Grid *grid = VectorGrid(vdummy_3d);
  int sg;

  ForSubgridI(sg, GridSubgrids(grid))
  {
    Subgrid *subgrid = GridSubgrid(grid, sg);

    int ix = SubgridIX(subgrid);
    int iy = SubgridIY(subgrid);
    int iz = SubgridIZ(subgrid);

    int nx = SubgridNX(subgrid);
    int ny = SubgridNY(subgrid);
    int nz = SubgridNZ(subgrid);

    for(i=0;i<nriverid;i++){
      if(riveridx[i]>=ix && riveridx[i]<(ix+nx) && riveridy[i]>=iy && riveridy[i]<(iy+ny)){
        for(j=iz;j<(iz+nz);j++){
          idx = nx*ny*j + (riveridy[i]-iy)*nx + (riveridx[i]-ix);
          if(pf_updateflag == 1) pf_statevec[idx] = subvec_p[idx];
          if(pf_updateflag == 2) pf_statevec[idx] = subvec_sat[idx]*subvec_porosity[idx];
          if(pf_updateflag == 3) pf_statevec[idx+enkf_subvecsize] = subvec_p[idx];
        }
      }
    }
  }
}

void init_n_domains_size(int* n_domains_p)
{
  int nshift = 0;
  /* state updates */
  // if(pf_updateflag == 1 || pf_updateflag == 2) {
    *n_domains_p = nx_local * ny_local;
  /*   nshift = nx_local * ny_local; */
  /* }   */
  /* else if(pf_updateflag == 3) { */
  /*   *n_domains_p = nx_local * ny_local * 2; */
  /*   nshift = 2 * nx_local * ny_local; */
  /* } */
  /* parameter updates   */
  /* if(pf_paramupdate == 1){ */
  /*   *n_domains_p = nshift + nx_local * ny_local; */
  /* } */
  /* else if(pf_paramupdate == 2){ */
  /*   *n_domains_p = nshift + pf_paramvecsize; */
  /* }  */
}

void init_parf_l_size(int* dim_l)
{
  int nshift = 0;
  /* state updates */
  if(pf_updateflag == 1) {
    *dim_l = nz_local;
    nshift = nz_local;
  }
  else if(pf_updateflag == 2) {
#ifdef FOR2131
    *dim_l = 2 * nz_local;
    nshift = 2 * nz_local;
#else
    *dim_l = nz_local;
    nshift = nz_local;
#endif
  }
  else if(pf_updateflag == 3) {
    *dim_l = 2 * nz_local;
    nshift = 2 * nz_local;
  }
  /* parameter updates */
  if(pf_paramupdate == 1){
    *dim_l = nshift + nz_local;
  }
  else if(pf_paramupdate == 2){
    *dim_l = nshift + 1;
  }
}
