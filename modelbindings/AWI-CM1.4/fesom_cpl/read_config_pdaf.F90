!$Id: read_config_pdaf.F90 2136 2019-11-22 18:56:35Z lnerger $
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
       ONLY: dim_state, dim_state_p, dim_ens, &
       offset, screen, filtertype, subtype, dim_ens, &
       delt_obs_ocn, delt_obs_atm, bias_obs, rms_obs, dim_bias, &
       incremental, type_forget, forget, &
       locweight, local_range, srange, loc_radius, &
       type_trans, type_sqrt, step_null, &
       path_obs_sst, file_sst_prefix, file_sst_suffix, &
       path_init, file_init, file_inistate, read_inistate, varscale, &
       sst_exclude_ice, sst_exclude_diff
  USE output_pdaf, &
       ONLY: write_da, write_ens, str_daspec

  IMPLICIT NONE
!EOP

! Local variables
  CHARACTER(len=100) :: nmlfile ='namelist.pdaf'    ! name of namelist file
  CHARACTER(len=32)  :: handle             ! Handle for command line parser
  LOGICAL :: printconfig = .TRUE.          ! Print information on all configuration parameters


  NAMELIST /pdaf/ filtertype, subtype, dim_ens, screen, &
       incremental, type_forget, forget, dim_bias, &
       local_range, locweight, srange, rms_obs, &
       path_obs_sst, file_sst_prefix, file_sst_suffix, &
       n_modeltasks, path_init, file_init, step_null, printconfig, &
       file_inistate, read_inistate, write_da, write_ens, varscale, &
       str_daspec, type_trans, type_sqrt, bias_obs, &
       delt_obs_ocn, delt_obs_atm, &     
       sst_exclude_ice, sst_exclude_diff
  

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
  CALL add_slash(path_init)

! *** Print configuration variables ***
  showconf: IF (printconfig .AND. mype_model==0 .AND. task_id==1) THEN

     WRITE (*,'(/a,1x,a)') 'FESOM-PDAF','-- Overview of PDAF configuration --'
     WRITE (*,'(a,3x,a)') 'FESOM-PDAF','PDAF [namelist: pdaf]:'
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','filtertype  ', filtertype
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','subtype     ', subtype
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','n_modeltasks', n_modeltasks
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','dim_ens     ', dim_ens
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','delt_obs_ocn', delt_obs_ocn
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','delt_obs_atm', delt_obs_atm
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','step_null   ', step_null
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','screen      ', screen
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','incremental ', incremental
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','type_forget ', type_forget
     WRITE (*,'(a,5x,a,f10.2)') 'FESOM-PDAF','forget      ', forget
     WRITE (*,'(a,5x,a,es10.2)') 'FESOM-PDAF','varscale    ', varscale
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','dim_bias    ', dim_bias
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','type_trans  ', type_trans
     WRITE (*,'(a,5x,a,es10.2)')'FESOM-PDAF','local_range ', local_range
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','locweight   ', locweight
     WRITE (*,'(a,5x,a,es10.2)') 'FESOM-PDAF','srange      ', srange
     WRITE (*,'(a,5x,a,es10.2)')'FESOM-PDAF','rms_obs     ', rms_obs
     WRITE (*,'(a,5x,a,es10.2)')'FESOM-PDAF','bias_obs    ', bias_obs
     WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','sst_exclude_ice', sst_exclude_ice
     WRITE (*,'(a,5x,a,f11.3)') 'FESOM-PDAF','sst_exclude_diff', sst_exclude_diff
     WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','path_obs_sst ', TRIM(path_obs_sst)
     WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','file_sst_prefix', TRIM(file_sst_prefix)
     WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','file_sst_suffix', TRIM(file_sst_suffix)
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
