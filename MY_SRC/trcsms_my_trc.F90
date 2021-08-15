MODULE trcsms_my_trc
   !!======================================================================
   !!                         ***  MODULE trcsms_my_trc  ***
   !! TOP :   Main module of the MY_TRC tracers
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) Original code
   !!----------------------------------------------------------------------
#if defined key_my_trc
   !!----------------------------------------------------------------------
   !!   'key_my_trc'                                               CFC tracers
   !!----------------------------------------------------------------------
   !! trc_sms_my_trc       : MY_TRC model main routine
   !! trc_sms_my_trc_alloc : allocate arrays specific to MY_TRC sms
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc         ! Ocean variables
   USE trc             ! TOP variables
   USE trd_oce
   USE trdtrc
   USE sbcflx_ano      ! Yona
   USE phycst          ! physical constants
   USE sbc_oce !,  ONLY : ln_flx_ano, nn_fsbc, nn_isf, fwfisf, rnf     ! Yona
   USE iom             ! IOM library !Yona
   USE sbcisf !Yona
   USE sbcrnf !Yona
   USE oce    !Yona
   USE traqsr, ONLY : nksr !Yona

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_my_trc       ! called by trcsms.F90 module
   PUBLIC   trc_sms_my_trc_alloc ! called by trcini_my_trc.F90 module

   ! Defined HERE the arrays specific to MY_TRC sms and ALLOCATE them in trc_sms_my_trc_alloc

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcsms_my_trc.F90 5385 2015-06-09 13:50:42Z cetlod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
      !!* Substitution
#  include "top_substitute.h90"
CONTAINS

   SUBROUTINE trc_sms_my_trc( kt )
      !!----------------------------------------------------------------------
      !!                     ***  trc_sms_my_trc  ***
      !!
      !! ** Purpose :   main routine of MY_TRC model
      !!
      !! ** Method  : -
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER ::  ji, jj, jk, jn   ! dummy loop index
      INTEGER ::  ii0, ii1, ij0, ij1   ! dummy loop index
      INTEGER ::  ikt, ikb
      REAL(wp)::  z1_e3t 
      REAL(wp), POINTER, DIMENSION(:,:,:) :: ztrmyt
      REAL(wp), DIMENSION(jpi,jpj) :: zhflxtot, zfwflxtot ! Heat and freshwater flux anomaly to force the tracers at the surface
      !and on the vertical for rnf, isf and qsr
      REAL(wp)                     ::   zt_frz            ! Freezing point temperature
      !!----------------------------------------------------------------------
      !i
      IF( nn_timing == 1 )  CALL timing_start('trc_sms_my_trc')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_sms_my_trc:  MY_TRC model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

      !----------------------------------------
      !       Surface fluxes
      !----------------------------------------

      zhflxtot(:,:) = 0.    ;    zfwflxtot(:,:) = 0.
      IF( ln_flx_ano ) THEN ! Always True
         IF( MOD( kt-1, nn_fsbc ) == 0 ) THEN                        ! update ocean fluxes at each SBC frequency
            ! Add sensible heat flux to heat heat flux from iceshelves
            zt_frz = -1.9 
            hfisf_ano(:,:) = hfisf_ano(:,:) - fwisf_ano(:,:) * zt_frz * rcp
         ENDIF

         zhflxtot(:,:) = qns_ano(:,:) ! Add qns at the surface to the trend; the other heat fluxes are distributed on vertical levels
         IF( ln_fwf_ano ) THEN ! freshwater fluxes anomalies already in hdiv
             ! Remove emp (from piC) from trend (emp=emp+emp_ano) to compensate for hdiv, to keep only effect of emp_ano 
            zfwflxtot(:,:) = sfx_ano(:,:) - (emp(:,:) - emp_ano(:,:)) * trb(:,:,1,jp_myt1)
         ELSE
             ! Remove emp from trend (but this time emp=emp in the physics because no fwf ano)
            zfwflxtot(:,:) = sfx_ano(:,:) + &
               &             (emp_ano(:,:)-emp(:,:)) * trb(:,:,1,jp_myt1) ! fwrnf_ano and fwisf_ano are distributed on vertical levels below
         ENDIF
      ENDIF
      ! Update tracer trends with flux anomalies at the surface only
      tra(:,:,1,jp_myt0) = tra(:,:,1,jp_myt0) + r1_rau0_rcp * zhflxtot(:,:) / fse3t(:,:,1) !ºC/s
      tra(:,:,1,jp_myt1) = tra(:,:,1,jp_myt1) + r1_rau0     * zfwflxtot(:,:) / fse3t(:,:,1) !0.001/s

      IF (ln_flx_ano) THEN
         
         !----------------------------------------
         !       Ice Shelf effects (ISF)
         ! Add iceshelf heat flux anomaly to the trend for PAT
         ! Add fwfisf anomaly in PAS trend only if it is not already read in the model physics (e.g. in CTL, STRESS and HEAT runs)
         !----------------------------------------
         !

         IF( nn_isf > 0 ) THEN
            DO jj = 1,jpj
               DO ji = 1,jpi
                  ikt = misfkt(ji,jj)
                  ikb = misfkb(ji,jj)
                  ! level fully included in the ice shelf boundary layer
                  DO jk = ikt, ikb - 1
                     tra(ji,jj,jk,jp_myt0) = tra(ji,jj,jk,jp_myt0)                                          &
                             &           + hfisf_ano(ji,jj) * r1_rau0_rcp * r1_hisf_tbl(ji,jj)
                     IF( .NOT. ln_fwf_ano ) THEN !Remove fwfisf from trend to compensate its presence in hdiv
                        tra(ji,jj,jk,jp_myt1) = tra(ji,jj,jk,jp_myt1) + &
                         &   (fwisf_ano(ji,jj)-fwfisf(ji,jj)) * trb(ji,jj,jk,jp_myt1) * r1_rau0 * r1_hisf_tbl(ji,jj)
                     ELSE !Idem but first remove fwisf_ano from fwfisf 
                        tra(ji,jj,jk,jp_myt1) = tra(ji,jj,jk,jp_myt1) - &
                         &   (fwfisf(ji,jj)-fwisf_ano(ji,jj)) * trb(ji,jj,jk,jp_myt1) * r1_rau0 * r1_hisf_tbl(ji,jj)
                     ENDIF
                  END DO
                  ! level partially included in ice shelf boundary layer
                  jk = ikb
                  tra(ji,jj,jk,jp_myt0) = tra(ji,jj,jk,jp_myt0)                                           &
                    &              + hfisf_ano(ji,jj) * r1_rau0_rcp * r1_hisf_tbl(ji,jj) * ralpha(ji,jj)
                  IF( .NOT. ln_fwf_ano ) THEN
                    tra(ji,jj,jk,jp_myt1) = tra(ji,jj,jk,jp_myt1) + &
                         &   (fwisf_ano(ji,jj)-fwfisf(ji,jj)) * trb(ji,jj,jk,jp_myt1) * r1_rau0 * r1_hisf_tbl(ji,jj) * ralpha(ji,jj)
                  ELSE  
                    tra(ji,jj,jk,jp_myt1) = tra(ji,jj,jk,jp_myt1) - &
                         &   (fwfisf(ji,jj)-fwisf_ano(ji,jj)) * trb(ji,jj,jk,jp_myt1) * r1_rau0 * r1_hisf_tbl(ji,jj) * ralpha(ji,jj)
                  ENDIF
               END DO
            END DO
         ENDIF

         !
         !----------------------------------------
         !        River Runoff effects
         ! Add runoffs heat flux anomaly to the trend for PAT 
         ! Add runoffs anomaly in PAS trend only if it is not already read in the model physics (e.g. in CTL, STRESS and HEAT runs)
         !----------------------------------------

         IF( ln_rnf ) THEN 
            DO jj = 1, jpj
               DO ji = 1, jpi
                  IF( fwrnf_ano(ji,jj) /= 0._wp ) THEN
                     DO jk = 1, nk_rnf(ji,jj)
                        tra(ji,jj,jk,jp_myt0) = tra(ji,jj,jk,jp_myt0)   &
                                         &   +  hfrnf_ano(ji,jj) * r1_rau0_rcp / h_rnf(ji,jj)
                        IF( .NOT. ln_fwf_ano ) THEN !Remove rnf from trend to compensate its presence in hdiv
                           tra(ji,jj,jk,jp_myt1) = tra(ji,jj,jk,jp_myt1) - &
                            &   (fwrnf_ano(ji,jj)-rnf(ji,jj)) * trb(ji,jj,jk,jp_myt1) * r1_rau0 / h_rnf(ji,jj)
                        ELSE !Idem but first remove fwrnf_ano from rnf
                           tra(ji,jj,jk,jp_myt1) = tra(ji,jj,jk,jp_myt1) + &
                            &   (rnf(ji,jj)-fwrnf_ano(ji,jj)) * trb(ji,jj,jk,jp_myt1) * r1_rau0 / h_rnf(ji,jj)
                        ENDIF
                     END DO
                  ENDIF
               END DO
            END DO
         ENDIF

         !----------------------------------------
         !        qsr
         ! Add qsr_ano to PAT trend distributed on the vertical as in traqsr
         !----------------------------------------
         DO jk = 1, nksr
           DO jj = 1, jpj
               DO ji = 1, jpi
                  z1_e3t = 1. / fse3t(ji,jj,jk)
                  tra(ji,jj,jk,jp_myt0) = tra(ji,jj,jk,jp_myt0) + qsr_ano(ji,jj) * fraqsr(ji,jj,jk) * r1_rau0_rcp * z1_e3t
               END DO
            END DO
         END DO

         !---------------------------------

         !---------------------------------------- 
         !Yona change trend values to 0 locally in a detroit in the Arctic where water column is very shallow and stratified
         !----------------------------------------

        ! ii0 = 208    ;   ii1 = 222 
        ! ij0 = 319    ;   ij1 = 329 
        ! tra( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1),:,jp_myt0 ) =  0. !PAT trend to 0
         !tra( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1),:,jp_myt1 ) =  0. !PAS trend to 0, then remove emp at the surface and rnf on the vertical

        ! !remove emp at the surface (negative sign to remove)
        ! IF( ln_fwf_ano ) THEN
        !     ! Remove emp (from piC) from trend (emp=emp+emp_ano) to compensate for hdiv, to keep only effect of emp_ano
        !    tra( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1),1,jp_myt1 ) = - (emp(mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1)) - &
        !        &   emp_ano(mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1))) * trb(mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1),1,jp_myt1) &
        !        &   * r1_rau0 / fse3t(mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1)
        ! ELSE
        !     ! Remove emp from trend (but this time emp=emp in the physics because no fwf ano)
        !     tra( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1),1,jp_myt1 ) = - emp(mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1)) * &
        !     &  trb(mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1,jp_myt1) * r1_rau0 / fse3t(mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),1)
        ! ENDIF
        ! !remove rnf on the vertical (positive sign to remove)
        ! IF( ln_rnf ) THEN
        !    DO jj = mj0(ij0), mj1(ij1)
        !       DO ji = mi0(ii0), mi1(ii1)
        !          IF( fwrnf_ano(ji,jj) /= 0._wp ) THEN
        !             DO jk = 1, nk_rnf(ji,jj)
        !                IF( .NOT. ln_fwf_ano ) THEN !Remove rnf from trend to compensate its presence in hdiv
        !                   tra(ji,jj,jk,jp_myt1) = tra(ji,jj,jk,jp_myt1) + &
        !                    &   rnf(ji,jj) * trb(ji,jj,jk,jp_myt1) * r1_rau0 / h_rnf(ji,jj)
        !                ELSE !Idem but first remove fwrnf_ano from rnf
        !                   tra(ji,jj,jk,jp_myt1) = tra(ji,jj,jk,jp_myt1) + &
        !                    &   (rnf(ji,jj)-fwrnf_ano(ji,jj)) * trb(ji,jj,jk,jp_myt1) * r1_rau0 / h_rnf(ji,jj)
        !                ENDIF
        !             END DO
        !          ENDIF
        !       END DO
        !    END DO
        ! ENDIF

         !----------------------------------------
         !

      ENDIF

!      CALL iom_put( "zhflxtot", zhflxtot(:,:) ) ! Output total heatflux anomaly seen by PAT trend
!      CALL iom_put( "zfwflxtot", zfwflxtot(:,:) ) ! Output total saltflux anomaly seen by PAS trend  
!      CALL iom_put( "tra_PAT", r1_rau0_rcp * zhflxtot(:,:) / fse3t(:,:,1) ) ! Output PAT trend
!      CALL iom_put( "tra_PAS", r1_rau0 * zfwflxtot(:,:) / fse3t(:,:,1) ) ! Output PAS trend

      IF( l_trdtrc )  CALL wrk_alloc( jpi, jpj, jpk, ztrmyt )

      IF( l_trdtrc ) THEN      ! Save the trends in the ixed layer
          DO jn = jp_myt0, jp_myt1
            ztrmyt(:,:,:) = tra(:,:,:,jn)
            CALL trd_trc( ztrmyt, jn, jptra_sms, kt )   ! save trends
          END DO
          CALL wrk_dealloc( jpi, jpj, jpk, ztrmyt )
      END IF
      !
      IF( nn_timing == 1 )  CALL timing_stop('trc_sms_my_trc')
      !
   END SUBROUTINE trc_sms_my_trc


   INTEGER FUNCTION trc_sms_my_trc_alloc()
      !!----------------------------------------------------------------------
      !!              ***  ROUTINE trc_sms_my_trc_alloc  ***
      !!----------------------------------------------------------------------
      !
      ! ALLOCATE here the arrays specific to MY_TRC
      ! ALLOCATE( tab(...) , STAT=trc_sms_my_trc_alloc )
      trc_sms_my_trc_alloc = 0      ! set to zero if no array to be allocated
      !
      IF( trc_sms_my_trc_alloc /= 0 ) CALL ctl_warn('trc_sms_my_trc_alloc : failed to allocate arrays')
      !
   END FUNCTION trc_sms_my_trc_alloc


#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                        No MY_TRC model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_sms_my_trc( kt )             ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'trc_sms_my_trc: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sms_my_trc
#endif

   !!======================================================================
END MODULE trcsms_my_trc
