!$Id$
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
       ONLY: mype_model, n_modeltasks, task_id, mype_filter_echam, &
       MPI_COMM_WORLD, MPI_INTEGER, MPIerr
  USE mod_assim_pdaf, &                          ! General variables for assimilation
       ONLY: dim_state, dim_state_p, dim_ens, dim_lag, dim_bias, &
       screen, step_null, filtertype, subtype, &
       DA_couple_type, incremental, type_trans, type_sqrt, &
       type_forget, forget, locweight, loctype, loc_ratio, &
       path_init, file_init, file_inistate, read_inistate, varscale, &
       use_global_obs, restart
  USE mod_assim_atm_pdaf, &                      ! General variables for assimilation into atmosphere
       ONLY: delt_obs_atm, delt_obs_atm_offset
  USE output_pdaf, &                             ! Variables for file output
       ONLY: write_da, write_ens, str_daspec
  USE obs_airt_pdafomi, &                        ! Variables for air temperature observations
       ONLY: assim_a_airt, rms_obs_airt, lradius_airt, sradius_airt

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

  ! Settings specific for the atmosphere
  NAMELIST /pdaf_atm/ screen, delt_obs_atm, delt_obs_atm_offset, &
       assim_a_airt, rms_obs_airt, lradius_airt, sradius_airt


! ****************************************************
! ***   Initialize PDAF parameters from namelist   ***
! ****************************************************

! *** Read namelist file ***
#ifdef DEBUG
  WRITE(*,*) 'Read PDAF namelist file: ',nmlfile
#endif

  OPEN (20,file=nmlfile)
  READ (20,NML=pdaf)
  READ (20,NML=pdaf_atm)
  CLOSE (20)

! *** Add trailing slash to paths ***

  CALL add_slash(path_init)

! *** Print configuration variables ***

  showconf: IF (printconfig .AND. mype_filter_echam==0 .AND. task_id==1) THEN

     WRITE (*,'(/a,1x,a)') 'ECHAM-PDAF','-- Overview of PDAF configuration --'
     WRITE (*,'(a,3x,a)') 'ECHAM-PDAF','PDAF [namelist: pdaf]:'
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','DA_couple_type  ', DA_couple_type
     WRITE (*,'(a,5x,a,l)')     'ECHAM-PDAF','restart     ', restart
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','filtertype  ', filtertype
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','subtype     ', subtype
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','n_modeltasks', n_modeltasks
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','dim_ens     ', dim_ens
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','delt_obs_atm_offset', delt_obs_atm_offset
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','delt_obs_atm', delt_obs_atm
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','step_null   ', step_null
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','screen      ', screen
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','incremental ', incremental
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','type_forget ', type_forget
     WRITE (*,'(a,5x,a,f10.2)') 'ECHAM-PDAF','forget      ', forget
     WRITE (*,'(a,5x,a,es10.2)') 'ECHAM-PDAF','varscale    ', varscale
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','dim_bias    ', dim_bias
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','type_trans  ', type_trans
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','locweight   ', locweight
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','loctype     ', loctype
     WRITE (*,'(a,5x,a,es10.2)') 'ECHAM-PDAF','loc_ratio   ', loc_ratio
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','dim_lag     ', dim_lag
     WRITE (*,'(a,5x,a,l)')     'ECHAM-PDAF','assim_a_airt ', assim_a_airt
     IF (assim_a_airt) THEN
        WRITE (*,'(a,5x,a,es10.2)')'ECHAM-PDAF','  rms_obs_airt ', rms_obs_airt
        WRITE (*,'(a,5x,a,es10.2)')'ECHAM-PDAF','  lradius_airt ', lradius_airt
        WRITE (*,'(a,5x,a,es10.2)')'ECHAM-PDAF','  sradius_airt', sradius_airt
     END IF
     WRITE (*,'(a,5x,a,a)')     'ECHAM-PDAF','path_init   ', TRIM(path_init)
     WRITE (*,'(a,5x,a,a)')     'ECHAM-PDAF','file_init   ', TRIM(file_init)
     IF (read_inistate) THEN
        WRITE (*,'(a,5x,a,a)')     'ECHAM-PDAF','file_inistate ', TRIM(file_inistate)
     ENDIF
     WRITE (*,'(a,5x,a,l)')     'ECHAM-PDAF','write_da    ', write_da
     WRITE (*,'(a,5x,a,l)')     'ECHAM-PDAF','write_ens   ', write_ens
     WRITE (*,'(a,5x,a,a)')     'ECHAM-PDAF','str_daspec  ',TRIM(str_daspec)
     WRITE (*,'(a,1x,a)') 'ECHAM-PDAF','-- End of PDAF configuration overview --'

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
