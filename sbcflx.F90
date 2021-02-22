MODULE sbcflx
   !!======================================================================
   !!                       ***  MODULE  sbcflx  ***
   !! Ocean forcing:  momentum, heat and freshwater flux formulation
   !!=====================================================================
   !! History :  1.0  !  2006-06  (G. Madec)  Original code
   !!            3.3  !  2010-10  (S. Masson)  add diurnal cycle
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   namflx   : flux formulation namlist
   !!   sbc_flx  : flux formulation as ocean surface boundary condition (forced mode, fluxes read in NetCDF files)
   !!----------------------------------------------------------------------
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
   USE sbcflx_ano !Yona

   IMPLICIT NONE
   PRIVATE

   PUBLIC sbc_flx       ! routine called by step.F90

   INTEGER , PARAMETER ::   jpfld   = 9   ! maximum number of files to read !Yona
   INTEGER , PARAMETER ::   jp_utau = 1   ! index of wind stress (i-component) file
   INTEGER , PARAMETER ::   jp_vtau = 2   ! index of wind stress (j-component) file
   INTEGER , PARAMETER ::   jp_qtot = 3   ! index of total (non solar+solar) heat file
   INTEGER , PARAMETER ::   jp_qsr  = 4   ! index of solar heat file
   INTEGER , PARAMETER ::   jp_emp  = 5   ! index of evaporation-precipation file
   INTEGER , PARAMETER ::   jp_sfx  = 6   ! index of salt flux flux !Yona
   INTEGER , PARAMETER ::   jp_sithic  = 7! index of sea ice thickness !Yona
   INTEGER , PARAMETER ::   jp_sivolu  = 8! index of sea ice volume per area !Yona
   INTEGER , PARAMETER ::   jp_siconc  = 9! index of sea ice fraction !Yona
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf    ! structure of input fields (file informations, fields read)

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO-consortium (2010)
   !! $Id: sbcflx.F90 4990 2014-12-15 16:42:49Z timgraham $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_flx( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_flx  ***
      !!
      !! ** Purpose :   provide at each time step the surface ocean fluxes
      !!                (momentum, heat, freshwater and runoff)
      !!
      !! ** Method  : - READ each fluxes in NetCDF files:
      !!                   i-component of the stress              utau  (N/m2)
      !!                   j-component of the stress              vtau  (N/m2)
      !!                   net downward heat flux                 qtot  (watt/m2)
      !!                   net downward radiative flux            qsr   (watt/m2)
      !!                   net upward freshwater (evapo - precip) emp   (kg/m2/s)
      !!                   salt flux                              sfx   (pss*h*rau/dt => g/m2/s) 
      !!                   sea ice thickness                      hcpl_i(m) !Yona
      !!                   sea ice volume per area                vcpl_i(m) !Yona
      !!                   sea ice fraction                       fr_i      !Yona
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
      !!              - hcpl_i      sea ice thickness !Yona
      !!              - vcpl_i      sea ice volume per area !Yona
      !!              - fr_i        sea ice fraction !Yona
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
      TYPE(FLD_N) ::   sn_utau, sn_vtau, sn_qtot, sn_qsr, sn_emp, sn_sfx, sn_sithic, sn_sivolu, sn_siconc  !Yona ! informations about the fields to be read
      NAMELIST/namsbc_flx/ cn_dir, sn_utau, sn_vtau, sn_qtot, sn_qsr, sn_emp, sn_sfx, sn_sithic, sn_sivolu, sn_siconc !Yona
      !!---------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN                ! First call kt=nit000
         ! set file information
         REWIND( numnam_ref )              ! Namelist namsbc_flx in reference namelist : Files for fluxes
         READ  ( numnam_ref, namsbc_flx, IOSTAT = ios, ERR = 901)
901      IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_flx in reference namelist', lwp )

         REWIND( numnam_cfg )              ! Namelist namsbc_flx in configuration namelist : Files for fluxes
         READ  ( numnam_cfg, namsbc_flx, IOSTAT = ios, ERR = 902 )
902      IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_flx in configuration namelist', lwp )
         IF(lwm) WRITE ( numond, namsbc_flx )
         !
         !                                         ! check: do we plan to use ln_dm2dc with non-daily forcing?
         IF( ln_dm2dc .AND. sn_qsr%nfreqh /= 24 )   &
            &   CALL ctl_stop( 'sbc_blk_core: ln_dm2dc can be activated only with daily short-wave forcing' )
         !
         !                                         ! store namelist information in an array
         slf_i(jp_utau) = sn_utau   ;   slf_i(jp_vtau) = sn_vtau !Yona
         slf_i(jp_qtot) = sn_qtot   ;   slf_i(jp_qsr ) = sn_qsr !Yona
         slf_i(jp_emp ) = sn_emp    ;   slf_i(jp_sfx ) = sn_sfx !Yona
         slf_i(jp_sithic) = sn_sithic !Yona
         slf_i(jp_sivolu) = sn_sivolu !Yona
         slf_i(jp_siconc) = sn_siconc !Yona

         !
         ALLOCATE( sf(jpfld), STAT=ierror )        ! set sf structure
         IF( ierror > 0 ) THEN
            CALL ctl_stop( 'sbc_flx: unable to allocate sf structure' )   ;   RETURN
         ENDIF
         DO ji= 1, jpfld
            ALLOCATE( sf(ji)%fnow(jpi,jpj,1) )
            IF( slf_i(ji)%ln_tint ) ALLOCATE( sf(ji)%fdta(jpi,jpj,1,2) )
         END DO
         !                                         ! fill sf with slf_i and control print
         CALL fld_fill( sf, slf_i, cn_dir, 'sbc_flx', 'flux formulation for ocean surface boundary condition', 'namsbc_flx' )
         !
      ENDIF

      CALL fld_read( kt, nn_fsbc, sf )                            ! input fields provided at the current time-step

      IF( MOD( kt-1, nn_fsbc ) == 0 ) THEN                        ! update ocean fluxes at each SBC frequency

         IF( ln_dm2dc ) THEN   ;   qsr(:,:) = sbc_dcy( sf(jp_qsr)%fnow(:,:,1) ) * tmask(:,:,1)   ! modify now Qsr to include the diurnal cycle
         ELSE                  ;   qsr(:,:) =          sf(jp_qsr)%fnow(:,:,1)   * tmask(:,:,1) !Yona add *tmask (reading fields from the coupled model)
         ENDIF
!CDIR COLLAPSE
         DO jj = 1, jpj                                           ! set the ocean fluxes from read fields
            DO ji = 1, jpi
               utau(ji,jj)   = sf(jp_utau)%fnow(ji,jj,1) * umask(ji,jj,1) !Yona mask fields coming from the coupled model
               vtau(ji,jj)   = sf(jp_vtau)%fnow(ji,jj,1) * vmask(ji,jj,1) !Yona
               qns (ji,jj)   = ( sf(jp_qtot)%fnow(ji,jj,1) - sf(jp_qsr)%fnow(ji,jj,1) ) * tmask(ji,jj,1) !Yona
               emp (ji,jj)   = sf(jp_emp )%fnow(ji,jj,1) * tmask(ji,jj,1) !Yona
               !sfx (ji,jj)  = sf(jp_sfx )%fnow(ji,jj,1) * 1000 *tmask(ji,jj,1) !Yona: convert input salflux (kg/m2/s) in g/m2/s
               sfx (ji,jj)   = sf(jp_sfx )%fnow(ji,jj,1) * tmask(ji,jj,1) !Yona: forcing already in correct unit g/m2/s
               hcpl_i(ji,jj) = sf(jp_sithic)%fnow(ji,jj,1) * tmask(ji,jj,1) !Yona
               vcpl_i(ji,jj) = sf(jp_sivolu)%fnow(ji,jj,1) * tmask(ji,jj,1) !Yona
               fr_i(ji,jj)   = sf(jp_siconc)%fnow(ji,jj,1) * tmask(ji,jj,1) !Yona
            END DO
         END DO
         !                                                        ! add to qns the heat due to e-p
         !!clem qns(:,:) = qns(:,:) - emp(:,:) * sst_m(:,:) * rcp        ! mass flux is at SST !Yona not used
         !
        CALL lbc_lnk_multi( utau, 'U', -1. , vtau, 'V', -1. ) !Yona: added lbc_lnk otherwise model blows up (bc fields from coupled model?)
        CALL lbc_lnk_multi( qns,'T',1., emp,'T',1., sfx,'T',1., qsr,'T',1., hcpl_i,'T',1., vcpl_i,'T',1., fr_i,'T',1. ) !Yona

         !                                                        ! module of wind stress and wind speed at T-point

         IF( nitend-nit000 <= 100 .AND. lwp ) THEN                ! control print (if less than 100 time-step asked)
            WRITE(numout,*)
            WRITE(numout,*) '        read daily momentum, heat and freshwater fluxes OK'
            DO jf = 1, jpfld
               IF( jf == jp_utau .OR. jf == jp_vtau )   zfact =     1.
               IF( jf == jp_qtot .OR. jf == jp_qsr  )   zfact =     0.1
               IF( jf == jp_emp                     )   zfact = 86400.
               WRITE(numout,*)
               WRITE(numout,*) ' day: ', ndastp , TRIM(sf(jf)%clvar), ' * ', zfact
               CALL prihre( sf(jf)%fnow, jpi, jpj, 1, jpi, 20, 1, jpj, 10, zfact, numout )
            END DO
            CALL FLUSH(numout)
         ENDIF
         !
      ENDIF

      !! Yona : add anomalies if they are activated
      IF( ln_flx_ano ) THEN
         CALL sbc_flx_ano( kt )
      ENDIF
      !! Yona
      !
      ! Note the use of 0.5*(2-umask) in order to unmask the stress along coastlines !!Yona not used anymore
      zcoef = 1. / ( zrhoa * zcdrag )
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vect. opt.
            ! ztx = ( utau(ji-1,jj  ) + utau(ji,jj) ) * 0.5 * ( 2. - MIN( umask(ji-1,jj  ,1), umask(ji,jj,1) ) )
            ! zty = ( vtau(ji  ,jj-1) + vtau(ji,jj) ) * 0.5 * ( 2. - MIN( vmask(ji  ,jj-1,1), vmask(ji,jj,1) ) )
            !!Yona, set like before otherwise it creates too many differences
            ztx = ( utau(ji-1,jj  ) + utau(ji,jj) )
            zty = ( vtau(ji  ,jj-1) + vtau(ji,jj) )
            !!Yona
            zmod = 0.5 * SQRT( ztx * ztx + zty * zty )
            taum(ji,jj) = zmod
            wndm(ji,jj) = SQRT( zmod * zcoef )  !!clem: not used
         END DO
      END DO
      taum(:,:) = taum(:,:) * tmask(:,:,1) ; wndm(:,:) = wndm(:,:) * tmask(:,:,1) !Yona
      CALL lbc_lnk( taum(:,:), 'T', 1. )   ;   CALL lbc_lnk( wndm(:,:), 'T', 1. ) !Yona
   END SUBROUTINE sbc_flx

   !!======================================================================
END MODULE sbcflx
