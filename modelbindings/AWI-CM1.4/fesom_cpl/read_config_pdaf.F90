!$Id: read_config_pdaf.F90 2395 2020-10-06 16:46:42Z lnerger $
!>  Routine to read configuration for PDAF for namelist
!!
!! This routine reads the namelist file with
!! parameters controlling data assimilation with PDAF.
!!
!! __Revision history:__
!! 2017-07 - Lars Nerger - Initial code for AWI-CM
!! * Later revisions - see repository log
!!
SUBROUTINE read_config_pdaf()

  USE mod_parallel_pdaf, &                       ! Variables for ensemble parallelization
       ONLY: mype_model, n_modeltasks, task_id, &
       MPI_COMM_WORLD, MPI_INTEGER, MPIerr
  USE mod_assim_pdaf, &                          ! General variables for assimilation
       ONLY: dim_state, dim_state_p, dim_ens, dim_lag, dim_bias, &
       screen, step_null, filtertype, subtype, &
       DA_couple_type, incremental, type_trans, type_sqrt, &
       type_forget, forget, locweight, loctype, loc_ratio, &
       path_init, file_init, file_inistate, read_inistate, varscale, &
       use_global_obs, restart
  USE mod_assim_oce_pdaf, &                      ! General variables for assimilation into ocean
       ONLY: delt_obs_ocn, delt_obs_ocn_offset
  USE output_pdaf, &                             ! Variables for file output
       ONLY: write_da, write_ens, str_daspec
  USE obs_SST_CMEMS_pdafomi, &                   ! Variables for CMEMS SST observations
       ONLY: assim_o_sst, path_obs_sst, file_sst_prefix, file_sst_suffix, &
       rms_obs_sst, bias_obs_sst, lradius_sst, sradius_sst, &
       sst_fixed_rmse, sst_exclude_diff, sst_exclude_ice, file_syntobs_sst 

  IMPLICIT NONE


! *** Local variables ***
  CHARACTER(len=100) :: nmlfile ='namelist.pdaf'    ! Name of namelist file
  CHARACTER(len=32)  :: handle                      ! Handle for command line parser
  LOGICAL :: printconfig = .TRUE.                   ! Print information on all configuration parameters


  ! General settings
  NAMELIST /pdaf/ n_modeltasks, dim_ens, dim_lag, dim_bias, filtertype, &
       subtype, incremental, type_forget, forget, &
       type_trans, type_sqrt, step_null, locweight, loctype, loc_ratio, &
       path_init, file_init, file_inistate, read_inistate, varscale, &
       write_da, write_ens, str_daspec, printconfig, &
       use_global_obs, restart

  ! Settings specific for the ocean
  NAMELIST /pdaf_oce/ screen, delt_obs_ocn, delt_obs_ocn_offset, &
       assim_o_sst, &                                                         ! SST
       path_obs_sst, file_sst_prefix, file_sst_suffix, &
       rms_obs_sst, bias_obs_sst, lradius_sst, sradius_sst,  &
       sst_exclude_ice, sst_exclude_diff, sst_fixed_rmse



! ****************************************************
! ***   Initialize PDAF parameters from namelist   ***
! ****************************************************

! *** Read namelist file ***
#ifdef DEBUG
  WRITE(*,*) 'Read PDAF namelist file: ',nmlfile
#endif

  OPEN (20,file=nmlfile)
  READ (20,NML=pdaf)
  READ (20,NML=pdaf_oce)
  CLOSE (20)

! *** Add trailing slash to paths ***
  CALL add_slash(path_obs_sst)
  CALL add_slash(path_init)

! *** Print configuration variables ***
  showconf: IF (printconfig .AND. mype_model==0 .AND. task_id==1) THEN

     WRITE (*,'(/a,1x,a)') 'FESOM-PDAF','-- Overview of PDAF configuration --'
     WRITE (*,'(a,3x,a)') 'FESOM-PDAF','PDAF [namelist: pdaf]:'
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','DA_couple_type  ', DA_couple_type
     WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','restart     ', restart
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','filtertype  ', filtertype
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','subtype     ', subtype
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','n_modeltasks', n_modeltasks
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','dim_ens     ', dim_ens
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','delt_obs_ocn_offset', delt_obs_ocn_offset
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','delt_obs_ocn', delt_obs_ocn
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','step_null   ', step_null
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','screen      ', screen
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','incremental ', incremental
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','type_forget ', type_forget
     WRITE (*,'(a,5x,a,f10.2)') 'FESOM-PDAF','forget      ', forget
     WRITE (*,'(a,5x,a,es10.2)') 'FESOM-PDAF','varscale    ', varscale
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','dim_bias    ', dim_bias
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','type_trans  ', type_trans
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','locweight   ', locweight
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','loctype     ', loctype
     WRITE (*,'(a,5x,a,es10.2)') 'FESOM-PDAF','loc_ratio   ', loc_ratio
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','dim_lag     ', dim_lag
     WRITE (*,'(a,5x,a,i8)')     'FESOM-PDAF','use_global_obs', use_global_obs
     WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','assim_o_sst   ', assim_o_sst
     IF (assim_o_sst) THEN
        WRITE (*,'(a,5x,a,es10.2)')'FESOM-PDAF','  rms_obs_sst ', rms_obs_sst
        WRITE (*,'(a,5x,a,es10.2)')'FESOM-PDAF','  bias_obs_sst  ', bias_obs_sst
        WRITE (*,'(a,5x,a,es10.2)')'FESOM-PDAF','  lradius_sst ', lradius_sst
        WRITE (*,'(a,5x,a,es10.2)')'FESOM-PDAF','  sradius_sst ', sradius_sst
        WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','  sst_fixed_rmse', sst_fixed_rmse
        WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','  sst_exclude_ice', sst_exclude_ice
        WRITE (*,'(a,5x,a,f11.3)') 'FESOM-PDAF','  sst_exclude_diff', sst_exclude_diff
        WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','  path_obs_sst     ', TRIM(path_obs_sst)
        WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','  file_sst_prefix  ', TRIM(file_sst_prefix)
        WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','  file_sst_suffix  ', TRIM(file_sst_suffix)
     END IF
     WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','path_init   ', TRIM(path_init)
     WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','file_init   ', TRIM(file_init)
     IF (read_inistate) THEN
        WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','file_inistate ', TRIM(file_inistate)
     ENDIF
     WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','write_da    ', write_da
     WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','write_ens   ', write_ens
     WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','str_daspec  ',TRIM(str_daspec)
     WRITE (*,'(a,1x,a)') 'FESOM-PDAF','-- End of PDAF configuration overview --'

  END IF showconf

END SUBROUTINE read_config_pdaf
! ==============================================================================
!BOP
!
! !ROUTINE: add_slash --- Add trailing slash to path string
!
! !INTERFACE:
SUBROUTINE add_slash(path)

! !DESCRIPTION:
! This routine ensures that a string defining a path
! has a trailing slash.
!
! !USES:
  IMPLICIT NONE

! !ARGUMENTS:
  CHARACTER(len=100) :: path  ! String holding the path
!EOP

! *** Local variables ***
  INTEGER :: strlength

! *** Add trailing slash ***
  strlength = LEN_TRIM(path)

  IF (path(strlength:strlength) /= '/') THEN
     path = TRIM(path) // '/'
  END IF

END SUBROUTINE add_slash
