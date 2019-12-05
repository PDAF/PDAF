!$Id: read_config_pdaf.f90 2135 2019-11-22 18:56:29Z lnerger $
!BOP
!
! !ROUTINE: read_config_pdaf - Read configuration for PDAF
!
! !INTERFACE: 
SUBROUTINE read_config_pdaf()

! !DESCRIPTION:
! This routine read the namelist file with
! parameters controlling data assimilation with
! PDAF.
!
! !REVISION HISTORY:
! 2017-07 - Lars Nerger - Initial code for AWI-CM
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel_pdaf, &
       ONLY: mype_model, n_modeltasks, task_id
  USE mod_assim_pdaf, & ! Variables for assimilation
       ONLY: dim_state, dim_state_p, dim_ens, dim_lag, &
       offset, screen, filtertype, subtype, dim_ens, &
       delt_obs_ocn, delt_obs_atm, &
       bias_obs, rms_obs, dim_bias, DA_couple_type, &
       incremental, type_forget, peak_obs_error, &
       forget, locweight, local_range, srange, &
       type_trans, type_sqrt, step_null, &
       eff_dim_obs, loc_radius, loctype, loc_ratio, &
       path_obs_sst, path_obs_prof, file_sst_prefix, &
       file_sst_suffix, file_prof_prefix, file_prof_suffix, &
       path_init, file_init, file_inistate, read_inistate, varscale, &
       rms_obs_T, rms_obs_S, writeprofile, &
       sst_exclude_ice, sst_exclude_diff, file_syntobs, twin_experiment, &
       dim_obs_max, prof_exclude_diff, &
       proffiles_o, assim_o_sst, assim_o_en4_t, assim_o_en4_s

  IMPLICIT NONE
!EOP

! Local variables
  CHARACTER(len=100) :: nmlfile ='namelist.pdaf'    ! name of namelist file
  CHARACTER(len=32)  :: handle             ! Handle for command line parser
  LOGICAL :: printconfig = .true.          ! Print information on all configuration parameters

  ! Temporary
  CHARACTER(len=100) :: str_daspec='DA'        ! String to identify assimilation experiment
  LOGICAL :: write_da = .true.                 ! Whether to write output file from assimilation
  LOGICAL :: write_ens = .false.               ! Whether to write ensemble files



  NAMELIST /pdaf/ filtertype, subtype, dim_ens, screen, &
       incremental, type_forget, forget, dim_bias, &
       local_range, locweight, srange, rms_obs, DA_couple_type, &
       path_obs_sst, path_obs_prof, file_sst_prefix, &
       file_sst_suffix, file_prof_prefix, file_prof_suffix, &
       n_modeltasks, peak_obs_error, &
       path_init, file_init, step_null, printconfig, &
       file_inistate, read_inistate, write_da, write_ens, varscale, &
       str_daspec, type_trans, type_sqrt, dim_lag, bias_obs, &
       loctype, loc_ratio, delt_obs_ocn, delt_obs_atm, &     
       rms_obs_T, rms_obs_S, writeprofile, &
       sst_exclude_ice, sst_exclude_diff, file_syntobs, twin_experiment, &
       dim_obs_max, prof_exclude_diff, proffiles_o, &
       assim_o_sst, assim_o_en4_t, assim_o_en4_s


! ****************************************************
! ***   Initialize PDAF parameters from namelist   ***
! ****************************************************

! *** Read namelist file ***
#ifdef DEBUG
  WRITE(*,*) 'Read PDAF namelist file: ',nmlfile
#endif

  OPEN (20,file=nmlfile)
  READ (20,NML=pdaf)
  CLOSE (20)

! *** Add trailing slash to paths ***
  CALL add_slash(path_obs_sst)
  CALL add_slash(path_obs_prof)
  CALL add_slash(path_init)

! *** Print configuration variables ***
  showconf: IF (printconfig .AND. mype_model==0 .AND. task_id==1) THEN

     WRITE (*,'(/a,1x,a)') 'ECHAM-PDAF','-- Overview of PDAF configuration --'
     WRITE (*,'(a,3x,a)') 'ECHAM-PDAF','PDAF [namelist: pdaf]:'
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','filtertype  ', filtertype
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','subtype     ', subtype
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','n_modeltasks', n_modeltasks
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','dim_ens     ', dim_ens
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','delt_obs_ocn', delt_obs_ocn
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','delt_obs_atm', delt_obs_atm
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','step_null   ', step_null
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','screen      ', screen
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','incremental ', incremental
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','type_forget ', type_forget
     WRITE (*,'(a,5x,a,f10.2)') 'ECHAM-PDAF','forget      ', forget
     WRITE (*,'(a,5x,a,es10.2)') 'ECHAM-PDAF','varscale    ', varscale
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','dim_bias    ', dim_bias
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','type_trans  ', type_trans
     WRITE (*,'(a,5x,a,es10.2)')'ECHAM-PDAF','local_range ', local_range
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','locweight   ', locweight
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','loctype     ', loctype
     WRITE (*,'(a,5x,a,es10.2)') 'ECHAM-PDAF','srange      ', srange
     WRITE (*,'(a,5x,a,es10.2)') 'ECHAM-PDAF','loc_ratio   ', loc_ratio
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','proffiles_o   ', proffiles_o
     WRITE (*,'(a,5x,a,l)')     'ECHAM-PDAF','assim_o_sst   ', assim_o_sst
     WRITE (*,'(a,5x,a,l)')     'ECHAM-PDAF','assim_o_en4_t ', assim_o_en4_t
     WRITE (*,'(a,5x,a,l)')     'ECHAM-PDAF','assim_o_en4_s ', assim_o_en4_s
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','writeprofile', writeprofile
     WRITE (*,'(a,5x,a,es10.2)')'ECHAM-PDAF','rms_obs     ', rms_obs
     WRITE (*,'(a,5x,a,es10.2)')'ECHAM-PDAF','rms_obs_T   ', rms_obs_T
     WRITE (*,'(a,5x,a,es10.2)')'ECHAM-PDAF','rms_obs_S   ', rms_obs_S
     WRITE (*,'(a,5x,a,es10.2)')'ECHAM-PDAF','peak_obs_error', peak_obs_error
     WRITE (*,'(a,5x,a,es10.2)')'ECHAM-PDAF','bias_obs    ', bias_obs
     WRITE (*,'(a,5x,a,l)')     'ECHAM-PDAF','sst_exclude_ice', sst_exclude_ice
     WRITE (*,'(a,5x,a,f11.3)') 'ECHAM-PDAF','sst_exclude_diff', sst_exclude_diff
     WRITE (*,'(a,5x,a,f11.3)') 'ECHAM-PDAF','prof_exclude_diff', prof_exclude_diff
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','dim_lag     ', dim_lag
     WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','DA_couple_type  ', DA_couple_type
     WRITE (*,'(a,5x,a,a)')     'ECHAM-PDAF','path_obs_sst ', TRIM(path_obs_sst)
     WRITE (*,'(a,5x,a,a)')     'ECHAM-PDAF','path_obs_prof', TRIM(path_obs_prof)
     WRITE (*,'(a,5x,a,a)')     'ECHAM-PDAF','file_sst_prefix', TRIM(file_sst_prefix)
     WRITE (*,'(a,5x,a,a)')     'ECHAM-PDAF','file_sst_suffix', TRIM(file_sst_suffix)
     WRITE (*,'(a,5x,a,a)')     'ECHAM-PDAF','file_prof_prefix', TRIM(file_prof_prefix)
     WRITE (*,'(a,5x,a,a)')     'ECHAM-PDAF','file_prof_suffix', TRIM(file_prof_suffix)
     WRITE (*,'(a,5x,a,a)')     'ECHAM-PDAF','path_init   ', TRIM(path_init)
     WRITE (*,'(a,5x,a,a)')     'ECHAM-PDAF','file_init   ', TRIM(file_init)
     IF (filtertype==11 .or. twin_experiment) THEN
        WRITE (*,'(a,5x,a,a)')     'ECHAM-PDAF','file_syntobs', TRIM(file_syntobs)
        WRITE (*,'(a,5x,a,i10)')   'ECHAM-PDAF','dim_obs_max ', dim_obs_max
     END IF
     IF (read_inistate) THEN
        WRITE (*,'(a,5x,a,a)')     'ECHAM-PDAF','file_inistate ', TRIM(file_inistate)
     ENDIF
     WRITE (*,'(a,5x,a,l)')     'ECHAM-PDAF','write_da    ', write_da
     WRITE (*,'(a,5x,a,l)')     'ECHAM-PDAF','write_ens   ', write_ens
     WRITE (*,'(a,5x,a,a)')     'ECHAM-PDAF','str_daspec  ',TRIM(str_daspec)
     WRITE (*,'(a,5x,a,l)')     'ECHAM-PDAF','twin_experiment', twin_experiment
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
