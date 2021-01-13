MODULE sbcflx_ano
   !!======================================================================
   !!                       ***  MODULE  sbcflx  ***
   !! Ocean forcing:  momentum, heat and freshwater flux formulation
   !!=====================================================================
   !! History :  1.0  !  2006-06  (G. Madec)  Original code
   !!            3.3  !  2010-10  (S. Masson)  add diurnal cycle
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   namflx_ano   : flux formulation namlist
   !!   sbc_flx_ano  : flux formulation as ocean surface boundary condition (forced mode, fluxes read in NetCDF files)
   !!----------------------------------------------------------------------
   !! Yona: adding this new routine to read the flux forcings perturbations 
   
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! surface boundary condition: ocean fields
   USE sbcdcy          ! surface boundary condition: diurnal cycle on qsr
   USE phycst          ! physical constants
   USE fldread         ! read input fields
   USE iom             ! IOM library
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   PUBLIC sbc_flx_ano       ! routine called by step.F90

   INTEGER , PARAMETER ::   jpfld   = 9  ! maximum number of files to read 
   INTEGER , PARAMETER ::   jp_utau = 1   ! index of wind stress (i-component) file
   INTEGER , PARAMETER ::   jp_vtau = 2   ! index of wind stress (j-component) file
   INTEGER , PARAMETER ::   jp_qtot = 3   ! index of total (non solar+solar) heat file
   INTEGER , PARAMETER ::   jp_qsr  = 4   ! index of solar heat file
   INTEGER , PARAMETER ::   jp_emp  = 5   ! index of evaporation-precipation file
   INTEGER , PARAMETER ::   jp_sfx  = 6   ! index of salt flux 
   INTEGER , PARAMETER ::   jp_hfrnf  = 7 ! index of rnf heat flux
   INTEGER , PARAMETER ::   jp_fwisf  = 8   ! index of iceshelf freshwater flux
   INTEGER , PARAMETER ::   jp_fwrnf  = 9   ! index of runoff heat flux
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf    ! structure of input fields (file informations, fields read)

   LOGICAL, PUBLIC     ::   ln_heat_ano        !: include heat flux anomalies or not (experiments ALL and HEAT)
   LOGICAL, PUBLIC     ::   ln_fwf_ano      !: include freshwater flux anomalies or not (experiments ALL and WATER)
   LOGICAL, PUBLIC     ::   ln_stress_ano      !: include wind stress anomalies or not (experiments ALL and STRESS)
   LOGICAL             ::   ln_dm2dc_ano      !: activate dirnal cycle for flux anomalies (monthly files) 

   INTEGER  :: nday_qsr_ano

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qsr_ano               !: solar heat flux anomaly
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qns_ano               !: non solar heat flux anomaly
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   emp_ano               !: emp anomaly
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   sfx_ano               !: salt flux anomaly under sea ice
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   utau_ano              !: zonal wind stress anomaly
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   vtau_ano              !: meridional wind stress anomaly
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   hfisf_ano             !: heat flux anomaly from iceshelf
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   fwisf_ano             !: freshwater flux anomaly from iceshelf
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   hfrnf_ano             !: heat flux anomlay of the runoffs
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   fwrnf_ano             !: freshwater flux anomaly from the runoffs

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO-consortium (2010) 
   !! $Id: sbcflx.F90 4990 2014-12-15 16:42:49Z timgraham $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_flx_ano( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_flx_ano  ***
      !!                   
      !! ** Purpose :   provide at each time step the surface ocean fluxes (Yona: perturbation terms)
      !!                (momentum, heat, freshwater and runoff) 
      !!
      !! ** Method  : - READ each fluxes in NetCDF files:
      !!                   i-component of the stress              utau_ano  (N/m2)
      !!                   j-component of the stress              vtau_ano  (N/m2)
      !!                   net downward heat flux                 qtot_ano  (watt/m2)
      !!                   net downward radiative flux            qsr_ano   (watt/m2)
      !!                   net upward freshwater (evapo - precip) emp_ano   (kg/m2/s)
      !!                   salt flux                              sfx_ano   (pss*h*rau/dt => g/m2/s)
      !!
      !!      CAUTION :  - never mask the surface stress fields
      !!                 - the stress is assumed to be in the (i,j) mesh referential
      !!
      !! ** Action  :   update at each time-step
      !!              - utau, vtau  i- and j-component of the wind stress
      !!              - taum        wind stress module at T-point
      !!              - wndm        10m wind module at T-point
      !!              - qns         non solar heat flux including heat flux due to emp
      !!              - qsr         solar heat flux
      !!              - emp         upward mass flux (empmr+runoffs)
      !!              - sfx         salt flux
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !!
      INTEGER  ::   ji, jj, jf            ! dummy indices
      INTEGER  ::   ierror                ! return error code
      INTEGER  ::   ios                   ! Local integer output status for namelist read
      REAL(wp) ::   zfact                 ! temporary scalar
      REAL(wp) ::   zrhoa  = 1.22         ! Air density kg/m3
      REAL(wp) ::   zcdrag = 1.5e-3       ! drag coefficient
      REAL(wp) ::   ztx, zty, zmod, zcoef ! temporary variables
      !!
      CHARACTER(len=100) ::  cn_dir                                       ! Root directory for location of flx files
      TYPE(FLD_N), DIMENSION(jpfld) ::   slf_i                            ! array of namelist information structures
      TYPE(FLD_N) ::   sn_utau, sn_vtau, sn_qtot, sn_qsr, sn_emp, sn_sfx,& 
              &  sn_hfrnf, sn_fwisf, sn_fwrnf  ! informations about the fields to be read
      NAMELIST/namsbc_flx_ano/ cn_dir, sn_utau, sn_vtau, sn_qtot, sn_qsr, sn_emp, &
              & sn_sfx, sn_hfrnf, sn_fwisf, sn_fwrnf, &   !Yona remettre bon commentaire
              & ln_heat_ano, ln_fwf_ano, ln_stress_ano, ln_dm2dc_ano
      !!---------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN                ! First call kt=nit000  
         ! set file information
         REWIND( numnam_ref )              ! Namelist namsbc_flx_ano in reference namelist : Files for fluxes
         READ  ( numnam_ref, namsbc_flx_ano, IOSTAT = ios, ERR = 901)
901      IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_flx_ano in reference namelist', lwp )

         REWIND( numnam_cfg )              ! Namelist namsbc_flx_ano in configuration namelist : Files for fluxes
         READ  ( numnam_cfg, namsbc_flx_ano, IOSTAT = ios, ERR = 902 )
902      IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_flx_ano in configuration namelist', lwp )
         IF(lwm) WRITE ( numond, namsbc_flx_ano ) 
         !
         !                                         ! check: do we plan to use ln_dm2dc with non-daily forcing?
         IF( ln_dm2dc_ano .AND. sn_qsr%nfreqh /= 24 )   &
            &   CALL ctl_stop( 'sbc_flx_ano: ln_dm2dc_ano can be activated only with daily short-wave forcing' ) 
         !
         IF( ln_dm2dc_ano )   nday_qsr_ano = -1   ! initialisation flag
         !                                         ! store namelist information in an array
         slf_i(jp_utau ) = sn_utau     ;   slf_i(jp_vtau ) = sn_vtau
         slf_i(jp_qtot ) = sn_qtot     ;   slf_i(jp_qsr  ) = sn_qsr 
         slf_i(jp_emp  ) = sn_emp      ;   slf_i(jp_sfx  ) = sn_sfx
         slf_i(jp_hfrnf) = sn_hfrnf
         slf_i(jp_fwisf) = sn_fwisf    ;   slf_i(jp_fwrnf) = sn_fwrnf
         !
         ALLOCATE( sf(jpfld), STAT=ierror )        ! set sf structure
         IF( ierror > 0 ) THEN   
            CALL ctl_stop( 'sbc_flx_ano: unable to allocate sf structure' )   ;   RETURN  
         ENDIF
         DO ji= 1, jpfld
            ALLOCATE( sf(ji)%fnow(jpi,jpj,1) )
            IF( slf_i(ji)%ln_tint ) ALLOCATE( sf(ji)%fdta(jpi,jpj,1,2) )
         END DO
         !                                         ! fill sf with slf_i and control print
         CALL fld_fill( sf, slf_i, cn_dir, 'sbc_flx_ano', 'flux formulation for ocean surface boundary condition (Yona:perturbation)', 'namsbc_flx_ano' )
         !
         ALLOCATE( qsr_ano(jpi,jpj)  , qns_ano(jpi,jpj)  , &
         &         emp_ano(jpi,jpj)  , sfx_ano(jpi,jpj)  , &
         &         utau_ano(jpi,jpj) , vtau_ano(jpi,jpj) , &
         &         hfisf_ano(jpi,jpj),  fwisf_ano(jpi,jpj), &
         &         hfrnf_ano(jpi,jpj),  fwrnf_ano(jpi,jpj), STAT=ierror )
         IF( ierror > 0 ) THEN   
            CALL ctl_stop( 'sbc_flx_ano: unable to allocate sf structure' )   ;   RETURN  
         ENDIF
         !
      ENDIF

      CALL fld_read( kt, nn_fsbc, sf )                            ! input fields provided at the current time-step
     
      IF( MOD( kt-1, nn_fsbc ) == 0 ) THEN                        ! update ocean fluxes at each SBC frequency
         IF( ln_dm2dc_ano ) THEN   ;   qsr_ano(:,:) = sbc_dcy( sf(jp_qsr)%fnow(:,:,1) * tmask(:,:,1), kday_qsr=nday_qsr_ano )   ! modify now Qsr to include the diurnal cycle
         ELSE                  ;   qsr_ano(:,:) =          sf(jp_qsr)%fnow(:,:,1) * tmask(:,:,1)
         ENDIF
!CDIR COLLAPSE
         DO jj = 1, jpj                                           ! set the ocean fluxes from read fields
            DO ji = 1, jpi
               utau_ano(ji,jj) = sf(jp_utau)%fnow(ji,jj,1) * umask(ji,jj,1)
               vtau_ano(ji,jj) = sf(jp_vtau)%fnow(ji,jj,1) * vmask(ji,jj,1)
               qns_ano (ji,jj) = (sf(jp_qtot)%fnow(ji,jj,1) - sf(jp_qsr)%fnow(ji,jj,1)) * tmask(ji,jj,1)
               emp_ano (ji,jj) = sf(jp_emp )%fnow(ji,jj,1) * tmask(ji,jj,1)
               sfx_ano (ji,jj) = sf(jp_sfx )%fnow(ji,jj,1) * tmask(ji,jj,1) !
               fwisf_ano (ji,jj) = -sf(jp_fwisf )%fnow(ji,jj,1) *  tmask(ji,jj,1)
               hfrnf_ano (ji,jj) = sf(jp_hfrnf )%fnow(ji,jj,1) * tmask(ji,jj,1)
               fwrnf_ano (ji,jj) = sf(jp_fwrnf )%fnow(ji,jj,1) * tmask(ji,jj,1)
            END DO
         END DO

         CALL lbc_lnk( utau_ano(:,:), 'U', -1. )
         CALL lbc_lnk( vtau_ano(:,:), 'V', -1. )
         CALL lbc_lnk( qns_ano (:,:), 'T',  1. )
         CALL lbc_lnk( qsr_ano (:,:), 'T',  1. )
         CALL lbc_lnk( emp_ano (:,:), 'T',  1. )
         CALL lbc_lnk( sfx_ano (:,:), 'T',  1. )
         CALL lbc_lnk( fwisf_ano (:,:), 'T',  1. )
         CALL lbc_lnk( hfrnf_ano (:,:), 'T',  1. )
         CALL lbc_lnk( fwrnf_ano (:,:), 'T',  1. )
         !                                                        ! add to qns the heat due to e-p
         !!clem qns(:,:) = qns(:,:) - emp(:,:) * sst_m(:,:) * rcp        ! mass flux is at SST
         !
         hfisf_ano(:,:) = fwisf_ano(:,:) * lfusisf !Yona compute hfisf from fwfisf

         !                                                 
         CALL iom_put( "emp_ano", emp_ano )                ! upward water flux
         CALL iom_put( "sfx_ano", sfx_ano  )                        ! downward salt flux  
         CALL iom_put( "qns_ano"   , qns_ano        )                   ! non solar heat flux
         CALL iom_put( "qsr_ano"   , qsr_ano  )                   ! solar heat flux
         CALL iom_put( "hfisf_ano" , hfisf_ano       )                   !  
         CALL iom_put( "hfrnf_ano"  , hfrnf_ano       )                   !
         CALL iom_put( "fwisf_ano"  , fwisf_ano       )                   !
         CALL iom_put( "fwrnf_ano"  , fwrnf_ano       )                   !
         !
      ENDIF
         CALL iom_put( "utau_ano", utau_ano )   ! i-wind stress   (stress can be updated at 
         CALL iom_put( "vtau_ano", vtau_ano )   ! j-wind stress    each time step in sea-ice)
      !
   END SUBROUTINE sbc_flx_ano

   !!======================================================================
END MODULE sbcflx_ano
