MODULE diahth
   !!======================================================================
   !!                       ***  MODULE  diahth  ***
   !! Ocean diagnostics: thermocline and 20 degree depth
   !!======================================================================
   !! History :  OPA  !  1994-09  (J.-P. Boulanger)  Original code
   !!                 !  1996-11  (E. Guilyardi)  OPA8 
   !!                 !  1997-08  (G. Madec)  optimization
   !!                 !  1999-07  (E. Guilyardi)  hd28 + heat content 
   !!            8.5  !  2002-06  (G. Madec)  F90: Free form and module
   !!   NEMO     3.2  !  2009-07  (S. Masson) hc300 bugfix + cleaning + add new diag
   !!   NEMO     3.6  !  2015_07  (J.F. Gueremy) hd26
   !!   NEMO     3.6  !  2016_03  (M. Chevallier ) hc700 + hc2000
   !!   NEMO     3.6  !  2017_09  (C. Ethe ) taking vvl into account properly + simplyfing + cleaning 
   !!----------------------------------------------------------------------
#if   defined key_diahth   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_diahth' :                              thermocline depth diag.
   !!----------------------------------------------------------------------
   !!   dia_hth      : Compute varius diagnostics associated with the thermal vertical structure
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library
   USE iom             ! I/O library
   USE timing          ! preformance summary

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dia_hth       ! routine called by step.F90
   PUBLIC   dia_hth_alloc ! routine called by nemogcm.F90

   LOGICAL , PUBLIC, PARAMETER          ::   lk_diahth = .TRUE.    !: thermocline-20d depths flag
   ! note: following variables should move to local variables once iom_put is always used 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hth    !: depth of the max vertical temperature gradient [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hd20   !: depth of 20 C isotherm                         [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hd26   !: depth of 26 C isotherm
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hd28   !: depth of 28 C isotherm                         [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   htc3   !: heat content of first 300 m                    [J/m2]
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   htc7   !: heat content of first 700 m                    [J/m2]
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   htc2   !: heat content of first 2000 m                   [J/m2]

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id$ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   FUNCTION dia_hth_alloc()
      !!---------------------------------------------------------------------
      INTEGER :: dia_hth_alloc
      !!---------------------------------------------------------------------
      !
      ALLOCATE(hth(jpi,jpj), hd20(jpi,jpj), hd26(jpi,jpj), hd28(jpi,jpj), htc3(jpi,jpj), &
                &  htc7(jpi,jpj),  htc2(jpi,jpj), STAT=dia_hth_alloc)
      !
      IF( lk_mpp           )   CALL mpp_sum ( dia_hth_alloc )
      IF(dia_hth_alloc /= 0)   CALL ctl_warn('dia_hth_alloc: failed to allocate arrays.')
      !
   END FUNCTION dia_hth_alloc


   SUBROUTINE dia_hth( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_hth  ***
      !!
      !! ** Purpose : Computes
      !!      the mixing layer depth (turbocline): avt = 5.e-4
      !!      the depth of strongest vertical temperature gradient
      !!      the mixed layer depth with density     criteria: rho = rho(10m or surf) + 0.03(or 0.01)
      !!      the mixed layer depth with temperature criteria: abs( tn - tn(10m) ) = 0.2       
      !!      the top of the thermochine: tn = tn(10m) - ztem2 
      !!      the pycnocline depth with density criteria equivalent to a temperature variation 
      !!                rho = rho10m + (dr/dT)(T,S,10m)*(-0.2 degC) 
      !!      the barrier layer thickness
      !!      the maximal verical inversion of temperature and its depth max( 0, max of tn - tn(10m) )
      !!      the depth of the 20 degree isotherm (linear interpolation)
      !!      the depth of the 26 degree isotherm (linear interpolation)
      !!      the depth of the 28 degree isotherm (linear interpolation)
      !!      the heat content of first 300 m (+ 700 m and 2000 m for CMIP6 :MCH)
      !!
      !! ** Method : 
      !!-------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      !!
      INTEGER                      ::   ji, jj, jk            ! dummy loop arguments
      REAL(wp)                     ::   zrho3 = 0.03_wp       ! density     criterion for mixed layer depth
      REAL(wp)                     ::   zrho1 = 0.01_wp       ! density     criterion for mixed layer depth
      REAL(wp)                     ::   ztem2 = 0.2_wp        ! temperature criterion for mixed layer depth
      REAL(wp)                     ::   zztmp, zzdep          ! temporary scalars inside do loop
      REAL(wp)                     ::   zu, zv, zw, zut, zvt  ! temporary workspace
      REAL(wp), DIMENSION(jpi,jpj) ::   zabs2      ! MLD: abs( tn - tn(10m) ) = ztem2 
      REAL(wp), DIMENSION(jpi,jpj) ::   ztm2       ! Top of thermocline: tn = tn(10m) - ztem2     
      REAL(wp), DIMENSION(jpi,jpj) ::   zrho10_3   ! MLD: rho = rho10m + zrho3      
      REAL(wp), DIMENSION(jpi,jpj) ::   zpycn      ! pycnocline: rho = rho10m + (dr/dT)(T,S,10m)*(-0.2 degC)
      REAL(wp), DIMENSION(jpi,jpj) ::   ztinv      ! max of temperature inversion
      REAL(wp), DIMENSION(jpi,jpj) ::   zdepinv    ! depth of temperature inversion
      REAL(wp), DIMENSION(jpi,jpj) ::   zrho0_3    ! MLD rho = rho(surf) = 0.03
      REAL(wp), DIMENSION(jpi,jpj) ::   zrho0_1    ! MLD rho = rho(surf) = 0.01
      REAL(wp), DIMENSION(jpi,jpj) ::   zmaxdzT    ! max of dT/dz
      REAL(wp), DIMENSION(jpi,jpj) ::   zdelr      ! delta rho equivalent to deltaT = 0.2
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   ztemp      ! temperature
      !!----------------------------------------------------------------------
      IF( nn_timing == 1 )   CALL timing_start('dia_hth')

      IF( kt == nit000 ) THEN
         !                                      ! allocate dia_hth array
         IF( dia_hth_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'lim_sbc_init : unable to allocate standard arrays' )
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dia_hth : diagnostics of the thermocline depth'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
         IF(lwp) WRITE(numout,*)
         !
      ENDIF

      ! initialization
      ztinv  (:,:) = 0._wp  
      zdepinv(:,:) = 0._wp  
      zmaxdzT(:,:) = 0._wp  
      DO jj = 1, jpj
         DO ji = 1, jpi
            zztmp = bathy(ji,jj)
            hth     (ji,jj) = zztmp
            zabs2   (ji,jj) = zztmp
            ztm2    (ji,jj) = zztmp
            zrho10_3(ji,jj) = zztmp
            zpycn   (ji,jj) = zztmp
        END DO
      END DO
      IF( nla10 > 1 ) THEN 
         DO jj = 1, jpj
            DO ji = 1, jpi
               zztmp = bathy(ji,jj)
               zrho0_3(ji,jj) = zztmp
               zrho0_1(ji,jj) = zztmp
            END DO
         END DO
      ENDIF
      
      ! Preliminary computation
      ! computation of zdelr = (dr/dT)(T,S,10m)*(-0.2 degC)
      DO jj = 1, jpj
         DO ji = 1, jpi
            IF( tmask(ji,jj,nla10) == 1. ) THEN
               zu  =  1779.50 + 11.250 * tsn(ji,jj,nla10,jp_tem) - 3.80   * tsn(ji,jj,nla10,jp_sal)                             &
                  &                                              - 0.0745 * tsn(ji,jj,nla10,jp_tem) * tsn(ji,jj,nla10,jp_tem)   &
                  &                                              - 0.0100 * tsn(ji,jj,nla10,jp_tem) * tsn(ji,jj,nla10,jp_sal)
               zv  =  5891.00 + 38.000 * tsn(ji,jj,nla10,jp_tem) + 3.00   * tsn(ji,jj,nla10,jp_sal)                             &
                  &                                              - 0.3750 * tsn(ji,jj,nla10,jp_tem) * tsn(ji,jj,nla10,jp_tem)
               zut =    11.25 -  0.149 * tsn(ji,jj,nla10,jp_tem) - 0.01   * tsn(ji,jj,nla10,jp_sal)
               zvt =    38.00 -  0.750 * tsn(ji,jj,nla10,jp_tem)
               zw  = (zu + 0.698*zv) * (zu + 0.698*zv)
               zdelr(ji,jj) = ztem2 * (1000.*(zut*zv - zvt*zu)/zw)
            ELSE
               zdelr(ji,jj) = 0._wp
            ENDIF
         END DO
      END DO

      ! ------------------------------------------------------------- !
      ! thermocline depth: strongest vertical gradient of temperature !
      ! turbocline depth (mixing layer depth): avt = zavt5            !
      ! MLD: rho = rho(1) + zrho3                                     !
      ! MLD: rho = rho(1) + zrho1                                     !
      ! ------------------------------------------------------------- !
      DO jk = jpkm1, 2, -1   ! loop from bottom to 2
         DO jj = 1, jpj
            DO ji = 1, jpi
               !
               zzdep = fsdepw(ji,jj,jk)
               zztmp = ( tsn(ji,jj,jk-1,jp_tem) - tsn(ji,jj,jk,jp_tem) ) / zzdep * tmask(ji,jj,jk)   ! vertical gradient of temperature (dT/dz)
               zzdep = zzdep * tmask(ji,jj,1)

               IF( zztmp > zmaxdzT(ji,jj) ) THEN                        
                  zmaxdzT(ji,jj) = zztmp   ;   hth    (ji,jj) = zzdep                ! max and depth of dT/dz
               ENDIF
               
               IF( nla10 > 1 ) THEN 
                  zztmp = rhop(ji,jj,jk) - rhop(ji,jj,1)                             ! delta rho(1)
                  IF( zztmp > zrho3 )          zrho0_3(ji,jj) = zzdep                ! > 0.03
                  IF( zztmp > zrho1 )          zrho0_1(ji,jj) = zzdep                ! > 0.01
               ENDIF

            END DO
         END DO
      END DO
      
      CALL iom_put( "mlddzt", hth )            ! depth of the thermocline
      IF( nla10 > 1 ) THEN 
         CALL iom_put( "mldr0_3", zrho0_3 )   ! MLD delta rho(surf) = 0.03
         CALL iom_put( "mldr0_1", zrho0_1 )   ! MLD delta rho(surf) = 0.01
      ENDIF

      ! ------------------------------------------------------------- !
      ! MLD: abs( tn - tn(10m) ) = ztem2                              !
      ! Top of thermocline: tn = tn(10m) - ztem2                      !
      ! MLD: rho = rho10m + zrho3                                     !
      ! pycnocline: rho = rho10m + (dr/dT)(T,S,10m)*(-0.2 degC)       !
      ! temperature inversion: max( 0, max of tn - tn(10m) )          !
      ! depth of temperature inversion                                !
      ! ------------------------------------------------------------- !
      DO jk = jpkm1, nlb10, -1   ! loop from bottom to nlb10
         DO jj = 1, jpj
            DO ji = 1, jpi
               !
               zzdep = fsdepw(ji,jj,jk) * tmask(ji,jj,1)
               !
               zztmp = tsn(ji,jj,nla10,jp_tem) - tsn(ji,jj,jk,jp_tem)  ! - delta T(10m)
               IF( ABS(zztmp) > ztem2 )      zabs2   (ji,jj) = zzdep   ! abs > 0.2
               IF(     zztmp  > ztem2 )      ztm2    (ji,jj) = zzdep   ! > 0.2
               zztmp = -zztmp                                          ! delta T(10m)
               IF( zztmp >  ztinv(ji,jj) ) THEN                        ! temperature inversion
                  ztinv(ji,jj) = zztmp   ;   zdepinv (ji,jj) = zzdep   ! max value and depth
               ENDIF

               zztmp = rhop(ji,jj,jk) - rhop(ji,jj,nla10)              ! delta rho(10m)
               IF( zztmp > zrho3        )    zrho10_3(ji,jj) = zzdep   ! > 0.03
               IF( zztmp > zdelr(ji,jj) )    zpycn   (ji,jj) = zzdep   ! > equi. delta T(10m) - 0.2
               !
            END DO
         END DO
      END DO

      CALL iom_put( "mld_dt02", zabs2        )   ! MLD abs(delta t) - 0.2
      CALL iom_put( "topthdep", ztm2         )   ! T(10) - 0.2
      CALL iom_put( "mldr10_3", zrho10_3     )   ! MLD delta rho(10m) = 0.03
      CALL iom_put( "pycndep" , zpycn        )   ! MLD delta rho equi. delta T(10m) = 0.2
      CALL iom_put( "tinv"    , ztinv        )   ! max. temp. inv. (t10 ref) 
      CALL iom_put( "depti"   , zdepinv      )   ! depth of max. temp. inv. (t10 ref) 


      ! ------------------------------- !
      !  Depth of 20C/26C/28C isotherm  !
      ! ------------------------------- !
      IF( iom_use ("20d") ) THEN  ! depth of the 20 isotherm
         ztem2 = 20.
         CALL dia_hth_dep( ztem2, hd20 )  
         CALL iom_put( "20d", hd20 )    
      ENDIF
      !
      IF( iom_use ("26d") ) THEN  ! depth of the 26 isotherm
         ztem2 = 26.
         CALL dia_hth_dep( ztem2, hd26 )  
         CALL iom_put( "26d", hd26 )    
      ENDIF
      !
      IF( iom_use ("28d") ) THEN  ! depth of the 28 isotherm
         ztem2 = 28.
         CALL dia_hth_dep( ztem2, hd28 )  
         CALL iom_put( "28d", hd28 )    
      ENDIF
        

      ! ----------------------------- !
      !  Heat content of first 300 m  !
      ! ----------------------------- !
      IF( iom_use ("hc300") .OR. iom_use ("hc700") .OR. iom_use ("hc2000") ) THEN  
         !ztemp(:,:,:) = tsn(:,:,:,jp_tem)  ! to get heat content in J/m2
         ztemp(:,:,:) = tsn(:,:,:,jp_tem) + rt0 ! to get heat content in K*m (CMIP6)
      ENDIF

      IF( iom_use ("hc300") ) THEN  
         zzdep = 300.
         CALL  dia_hth_htc( zzdep, ztemp, htc3 )
         !htc3(:,:) = rau0_rcp * htc3(:,:)  ! to get heat content in J/m2
         CALL iom_put( "hc300", htc3 )  ! vertically integrated heat content (in K*m)
      ENDIF
      !
      ! ----------------------------- !
      !  Heat content of first 700 m  !
      ! ----------------------------- !
      IF( iom_use ("hc700") ) THEN  
         zzdep = 700.
         CALL  dia_hth_htc( zzdep, ztemp, htc7 )
         !htc7(:,:) = rau0_rcp * htc7(:,:)  ! to get heat content in J/m2
         CALL iom_put( "hc700", htc7 )  ! vertically integrated heat content (K*m) 
      ENDIF
      !
      ! ----------------------------- !
      !  Heat content of first 2000 m  !
      ! ----------------------------- !
      IF( iom_use ("hc2000") ) THEN  
         zzdep = 2000.
         CALL  dia_hth_htc( zzdep, ztemp, htc2 )
         !htc2(:,:) = rau0_rcp * htc2(:,:)  ! to get heat content in J/m2
         CALL iom_put( "hc2000", htc2 )  ! vertically integrated heat content (K*m) 
      ENDIF
      !
      IF( nn_timing == 1 )   CALL timing_stop('dia_hth')
      !
   END SUBROUTINE dia_hth
   
   SUBROUTINE dia_hth_dep( ptem, pdept )
      !
      REAL(wp), INTENT(in) :: ptem
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pdept     
      !
      INTEGER  :: ji, jj, jk, iid
      REAL(wp) :: zztmp, zzdep
      INTEGER, DIMENSION(jpi,jpj) :: iktem
      
      ! --------------------------------------- !
      ! search deepest level above ptem         !
      ! --------------------------------------- !
      iktem(:,:) = 1
      DO jk = 1, jpkm1   ! beware temperature is not always decreasing with depth => loop from top to bottom
         DO jj = 1, jpj
            DO ji = 1, jpi
               zztmp = tsn(ji,jj,jk,jp_tem)
               IF( zztmp >= ptem )   iktem(ji,jj) = jk
            END DO
         END DO
      END DO

      ! ------------------------------- !
      !  Depth of ptem isotherm         !
      ! ------------------------------- !
      DO jj = 1, jpj
         DO ji = 1, jpi
            !
            zzdep = fsdepw(ji,jj,mbkt(ji,jj)+1)       ! depth of the ocean bottom
            !
            iid = iktem(ji,jj)
            IF( iid /= 1 ) THEN 
                zztmp =     fsdept(ji,jj,iid  )   &                     ! linear interpolation
                  &  + (    fsdept(ji,jj,iid+1) - fsdept(ji,jj,iid)                       )   &
                  &  * ( 20.*tmask(ji,jj,iid+1) - tsn(ji,jj,iid,jp_tem)                       )   &
                  &  / ( tsn(ji,jj,iid+1,jp_tem) - tsn(ji,jj,iid,jp_tem) + (1.-tmask(ji,jj,1)) )
               pdept(ji,jj) = MIN( zztmp , zzdep) * tmask(ji,jj,1)       ! bound by the ocean depth
            ELSE 
               pdept(ji,jj) = 0._wp
            ENDIF
         END DO
      END DO
      !
   END SUBROUTINE dia_hth_dep


   SUBROUTINE dia_hth_htc( pdep, ptn, phtc )
      !
      REAL(wp), INTENT(in) :: pdep     ! depth over the heat content
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in) :: ptn   
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: phtc  
      !
      INTEGER  :: ji, jj, jk, ik
      REAL(wp), DIMENSION(jpi,jpj) :: zthick
      INTEGER , DIMENSION(jpi,jpj) :: ilevel


      ! surface boundary condition
      IF( lk_vvl ) THEN   ;   zthick(:,:) = 0._wp       ;   phtc(:,:) = 0._wp                                   
      ELSE                ;   zthick(:,:) = sshn(:,:)   ;   phtc(:,:) = ptn(:,:,1) * sshn(:,:) * tmask(:,:,1)   
      ENDIF
      !
      ilevel(:,:) = 1
      DO jk = 2, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( ( fsdept(ji,jj,jk ) < pdep ) .AND. ( tmask(ji,jj,jk) == 1 ) ) THEN
                   ilevel(ji,jj) = jk
                   zthick(ji,jj) = zthick(ji,jj) + fse3t(ji,jj,jk)
                   phtc  (ji,jj) = phtc  (ji,jj) + fse3t(ji,jj,jk) * ptn(ji,jj,jk)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            ik = ilevel(ji,jj)
            zthick(ji,jj) = pdep - zthick(ji,jj)   !   remaining thickness to reach depht pdep
            phtc(ji,jj)   = phtc(ji,jj) + ptn(ji,jj,ik+1) * MIN( fse3t(ji,jj,ik+1), zthick(ji,jj) ) &
                                                                 * tmask(ji,jj,ik+1)
         END DO
      ENDDO
      !
      !
   END SUBROUTINE dia_hth_htc
   

#else
   !!----------------------------------------------------------------------
   !!   Default option :                                       Empty module
   !!----------------------------------------------------------------------
   LOGICAL , PUBLIC, PARAMETER ::   lk_diahth = .FALSE.  !: thermocline-20d depths flag
CONTAINS
   SUBROUTINE dia_hth( kt )         ! Empty routine
      WRITE(*,*) 'dia_hth: You should not have seen this print! error?', kt
   END SUBROUTINE dia_hth
#endif

   !!======================================================================
END MODULE diahth
