MODULE diaar5
   !!======================================================================
   !!                       ***  MODULE  diaar5  ***
   !! AR5 diagnostics
   !!======================================================================
   !! History :  3.2  !  2009-11  (S. Masson)  Original code
   !!            3.3  !  2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase + merge TRC-TRA
   !!----------------------------------------------------------------------
#if defined key_diaar5   || defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_diaar5'  :                           activate ar5 diagnotics
   !!----------------------------------------------------------------------
   !!   dia_ar5       : AR5 diagnostics
   !!   dia_ar5_init  : initialisation of AR5 diagnostics
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and active tracers 
   USE dom_oce        ! ocean space and time domain
   USE eosbn2         ! equation of state                (eos_bn2 routine)
   USE lib_mpp        ! distribued memory computing library
   USE lib_fortran     ! Fortran routines library 
   USE iom            ! I/O manager library
   USE timing         ! preformance summary
   USE wrk_nemo       ! working arrays
   USE fldread        ! type FLD_N
   USE phycst         ! physical constant
   USE in_out_manager  ! I/O manager
   USE zdfddm
   USE zdf_oce

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dia_ar5        ! routine called in step.F90 module
   PUBLIC   dia_ar5_init   ! routine called in opa.F90 module
   PUBLIC   dia_ar5_alloc  ! routine called in nemogcm.F90 module

   LOGICAL, PUBLIC, PARAMETER :: lk_diaar5 = .TRUE.   ! coupled flag

   REAL(wp)                         ::   vol0      ! ocean volume (interior domain)
   REAL(wp)                         ::   area_tot     ! total ocean surface (interior domain)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:  ) ::   area         ! cell surface (interior domain)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:  ) ::   thick0       ! ocean thickness (interior domain)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sn0          ! initial salinity
      
   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "zdfddm_substitute.h90"

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   FUNCTION dia_ar5_alloc()
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dia_ar5_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER :: dia_ar5_alloc
      !!----------------------------------------------------------------------
      !
      ALLOCATE( area(jpi,jpj), thick0(jpi,jpj), sn0(jpi,jpj,jpk) , STAT=dia_ar5_alloc )
      !
      IF( lk_mpp             )   CALL mpp_sum ( dia_ar5_alloc )
      IF( dia_ar5_alloc /= 0 )   CALL ctl_warn('dia_ar5_alloc: failed to allocate arrays')
      !
   END FUNCTION dia_ar5_alloc


   SUBROUTINE dia_ar5( kt )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dia_ar5  ***
      !!
      !! ** Purpose :   compute and output some AR5 diagnostics
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !
      INTEGER  ::   ji, jj, jk, iks, ikb                      ! dummy loop arguments
      REAL(wp) ::   zvolssh, zvol, zssh_steric, zztmp, zarho, ztemp, zsal, zmass
      REAL(wp) ::   zaw, zbw, zrw, zsst, zsss
      !
      REAL(wp), POINTER, DIMENSION(:,:)     :: zarea_ssh, zbotpres, zpe, z2d       ! 2D workspace 
      REAL(wp), POINTER, DIMENSION(:,:,:)   :: zrhd , zrhop, z3d, ztpot            ! 3D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:,:) :: ztsn                       ! 4D workspace
      !!--------------------------------------------------------------------
      IF( nn_timing == 1 )   CALL timing_start('dia_ar5')
 
      CALL wrk_alloc( jpi, jpj           , zarea_ssh, zbotpres, zpe, z2d )
      CALL wrk_alloc( jpi, jpj, jpk      , zrhd     , zrhop, z3d, ztpot  )
      CALL wrk_alloc( jpi, jpj, jpk, jpts, ztsn                          )

      !Call to init moved to here so that we can call iom_use in the
      !initialisation
      IF( kt == nit000 )     CALL dia_ar5_init

      zarea_ssh(:,:) = area(:,:) * sshn(:,:)
      !                                         ! total volume of liquid seawater
      zvolssh = glob_sum( zarea_ssh(:,:) ) 
      zvol    = vol0 + zvolssh
      !
      z3d(:,:,jpk) = 0._wp        ! ocean volume
      DO jk = 1, jpkm1
         z3d (:,:,jk) = area(:,:) * fse3t(:,:,jk) * tmask(:,:,jk)
      END DO

      CALL iom_put( "bathy"     , bathy(:,:) )
      CALL iom_put( "e2u", e2u(:,:) )
      CALL iom_put( "e1v", e1v(:,:) )
      CALL iom_put( "areacello" , area(:,:) )
      !
      CALL iom_put( "volcello"  , z3d(:,:,:)         )     ! WARNING not consistent with CMIP DR where volcello is at ca. 2000 
      CALL iom_put( "masscello" , rau0 * fse3t(:,:,:) * tmask(:,:,:) )  ! ocean mass 
      !
      CALL iom_put( 'voltot', zvol               )
      CALL iom_put( 'sshtot', zvolssh / area_tot )
      CALL iom_put( 'sshdyn', sshn(:,:) - (zvolssh / area_tot) )
      
      !                     
      IF( iom_use('sshthster')) THEN
         ztsn(:,:,:,jp_tem) = tsn(:,:,:,jp_tem)                    ! thermosteric ssh
         ztsn(:,:,:,jp_sal) = sn0(:,:,:)
         CALL eos( ztsn, zrhd, fsdept_n(:,:,:) )                       ! now in situ density using initial salinity
         !
         zbotpres(:,:) = 0._wp                        ! no atmospheric surface pressure, levitating sea-ice
         DO jk = 1, jpkm1
            zbotpres(:,:) = zbotpres(:,:) + fse3t(:,:,jk) * zrhd(:,:,jk)
         END DO
         IF( .NOT.lk_vvl ) THEN
            IF ( ln_isfcav ) THEN
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     iks = mikt(ji,jj)
                     zbotpres(ji,jj) = zbotpres(ji,jj) + sshn(ji,jj) * zrhd(ji,jj,iks) + riceload(ji,jj)
                  END DO
               END DO
            ELSE
              zbotpres(:,:) = zbotpres(:,:) + sshn(:,:) * zrhd(:,:,1)
           END IF
         END IF
         !                                         
         zarho = glob_sum( area(:,:) * zbotpres(:,:) ) 
         zssh_steric = - zarho / area_tot
         CALL iom_put( 'sshthster', zssh_steric )
        !
      ENDIF    
      !                                         ! steric sea surface height
      CALL eos( tsn, zrhd, zrhop, fsdept_n(:,:,:) )                 ! now in situ and potential density
      zbotpres(:,:) = 0._wp                        ! no atmospheric surface pressure, levitating sea-ice
      DO jk = 1, jpkm1
         zbotpres(:,:) = zbotpres(:,:) + fse3t(:,:,jk) * zrhd(:,:,jk)
      END DO
      IF( .NOT.lk_vvl ) THEN
         IF ( ln_isfcav ) THEN
            DO jj = 1, jpj
               DO ji = 1, jpi
                  iks = mikt(ji,jj)
                  zbotpres(ji,jj) = zbotpres(ji,jj) + sshn(ji,jj) * zrhd(ji,jj,iks) + riceload(ji,jj)
               END DO
            END DO
         ELSE
            zbotpres(:,:) = zbotpres(:,:) + sshn(:,:) * zrhd(:,:,1)
         END IF
      END IF
      !    
      zrhop(:,:,jpk) = 0._wp
      CALL iom_put( 'rhop', zrhop )
      !
      zarho = glob_sum( area(:,:) * zbotpres(:,:) ) 
      zssh_steric = - zarho / area_tot
      CALL iom_put( 'sshsteric', zssh_steric )
      
      !                                         ! ocean bottom pressure
      zztmp = rau0 * grav * 1.e-4_wp               ! recover pressure from pressure anomaly and cover to dbar = 1.e4 Pa
      zbotpres(:,:) = zztmp * ( zbotpres(:,:) + sshn(:,:) + thick0(:,:) )
      CALL iom_put( 'botpres', zbotpres )

      !                                         ! Mean density anomalie, temperature and salinity
      !
      ztsn(:,:,:,:) = 0._wp                    ! ztsn(:,:,1,jp_tem/sal) is used here as 2D Workspace for temperature & salinity
      DO jk = 1, jpkm1
         ztsn(:,:,1,jp_tem) = ztsn(:,:,1,jp_tem) + area(:,:) * fse3t(:,:,jk) * tsn(:,:,jk,jp_tem)
         ztsn(:,:,1,jp_sal) = ztsn(:,:,1,jp_sal) + area(:,:) * fse3t(:,:,jk) * tsn(:,:,jk,jp_sal)
      ENDDO
      IF( .NOT.lk_vvl ) THEN
         IF( ln_isfcav ) THEN
            DO jj = 1, jpj
               DO ji = 1, jpi
                  iks = mikt(ji,jj)
                  ztsn(ji,jj,1,jp_tem) = ztsn(ji,jj,1,jp_tem) + zarea_ssh(ji,jj) * tsn(ji,jj,iks,jp_tem) 
                  ztsn(ji,jj,1,jp_sal) = ztsn(ji,jj,1,jp_sal) + zarea_ssh(ji,jj) * tsn(ji,jj,iks,jp_sal) 
               END DO
            END DO
         ELSE
            ztsn(:,:,1,jp_tem) = ztsn(:,:,1,jp_tem) + zarea_ssh(:,:) * tsn(:,:,1,jp_tem) 
            ztsn(:,:,1,jp_sal) = ztsn(:,:,1,jp_sal) + zarea_ssh(:,:) * tsn(:,:,1,jp_sal) 
         END IF
      ENDIF
      !
      zsss  = glob_sum( area(:,:) * tsn(:,:,1,jp_sal) )
      ztemp = glob_sum( ztsn(:,:,1,jp_tem) )
      zsal  = glob_sum( ztsn(:,:,1,jp_sal) )
      zmass = rau0 * ( zarho + zvol )      
      !
      CALL iom_put( 'masstot', zmass             )          ! total mass of liquid seawater
      CALL iom_put( 'temptot', ztemp / zvol      )          ! potential temperature in liquid seawater
      CALL iom_put( 'saltot' , zsal  / zvol      )          ! Salinity of liquid seawater
      CALL iom_put( 'saltctot',zsal * rau0 *0.001)          ! Salt content of liquid seawater in kg of salt 
      CALL iom_put( 'ssstot' , zsss  / area_tot  )          ! Salinity of liquid seawater at surface

      IF( iom_use( "e3tb" ) )  THEN
         DO jj = 1, jpj
            DO ji = 1, jpi
               ikb = mbkt(ji,jj)
               z2d(ji,jj) = fse3t(ji,jj,ikb)
            END DO
         END DO
         CALL iom_put( "e3tb", z2d )  
      ENDIF 


 

      IF( nn_eos == -1 ) THEN        ! ! potential temperature (TEOS-10 case)
         ztpot(:,:,jpk) = 0._wp
         ztpot(:,:,:)   = eos_pt_from_ct( tsn(:,:,:,jp_tem), tsn(:,:,:,jp_sal) )
         !
         CALL iom_put( "toce_pot", ztpot(:,:,:) )  ! potential temperature (TEOS-10 case)
         CALL iom_put( "sst_pot" , ztpot(:,:,1) )  ! surface temperature
         !
         IF( iom_use('temptot_pot') ) THEN   ! Output potential temperature in case we use TEOS-10
            z2d(:,:) = 0._wp
            DO jk = 1, jpkm1
               z2d(:,:) = z2d(:,:) + area(:,:) * fse3t(:,:,jk) * ztpot(:,:,jk)
            END DO
            zsst  = glob_sum( area(:,:) * ztpot(:,:,1)  ) 
            ztemp = glob_sum( z2d(:,:)  ) 
            CALL iom_put( 'temptot_pot', ztemp / zvol )
            CALL iom_put( 'ssttot'     , zsst / area_tot )
         ENDIF
         ! Vertical integral of temperature
         IF( iom_use("tosmint_pot") ) THEN
           z2d(:,:) = 0._wp
           DO jk = 1, jpkm1
              DO jj = 1, jpj
                 DO ji = 1, jpi   ! vector opt.
                    z2d(ji,jj) = z2d(ji,jj) + rau0 * fse3t(ji,jj,jk) *  ztpot(ji,jj,jk)
                 END DO
              END DO
           END DO
           CALL iom_put( "tosmint_pot", z2d ) 
        ENDIF
      ELSE       
         zsst  = glob_sum( area(:,:) * tsn(:,:,1,jp_tem) ) ! Case EOS-80 : compute sst anyway
         CALL iom_put('ssttot', zsst / area_tot )
      ENDIF

      IF( iom_use( 'tnpeo' )) THEN    
      ! Work done against stratification by vertical mixing
      ! Exclude points where rn2 is negative as convection kicks in here and
      ! work is not being done against stratification
         zpe(:,:) = 0._wp
         IF( lk_zdfddm ) THEN
            DO jk = 2, jpk
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     IF( rn2(ji,jj,jk) > 0._wp ) THEN
                        zrw =   ( fsdepw(ji,jj,jk  ) - fsdept(ji,jj,jk) )   &
                           &  / ( fsdept(ji,jj,jk-1) - fsdept(ji,jj,jk) )
                        !
                        zaw = rab_n(ji,jj,jk,jp_tem) * (1. - zrw) + rab_n(ji,jj,jk-1,jp_tem)* zrw
                        zbw = rab_n(ji,jj,jk,jp_sal) * (1. - zrw) + rab_n(ji,jj,jk-1,jp_sal)* zrw
                        !
                        zpe(ji, jj) = zpe(ji, jj)            &
                           &        -  grav * (    avt(ji,jj,jk) * zaw * (tsn(ji,jj,jk-1,jp_tem) - tsn(ji,jj,jk,jp_tem) )  &
                           &                   - fsavs(ji,jj,jk) * zbw * (tsn(ji,jj,jk-1,jp_sal) - tsn(ji,jj,jk,jp_sal) ) )
                     ENDIF
                  END DO
               END DO
             END DO
          ELSE
            DO jk = 1, jpk
               DO ji = 1, jpi
                  DO jj = 1, jpj
                     zpe(ji,jj) = zpe(ji,jj) + avt(ji, jj, jk) * MIN(0._wp,rn2(ji, jj, jk)) * rau0 * fse3w(ji, jj, jk)
                  END DO
               END DO
            END DO
         ENDIF
         CALL lbc_lnk(zpe, 'T', 1._wp)         
         CALL iom_put( 'tnpeo', zpe  )
      ENDIF
       !
      CALL wrk_dealloc( jpi, jpj           , zarea_ssh, zbotpres, zpe, z2d )
      CALL wrk_dealloc( jpi, jpj, jpk      , zrhd     , zrhop, z3d, ztpot  )
      CALL wrk_dealloc( jpi, jpj, jpk, jpts, ztsn                          )
      !
      IF( nn_timing == 1 )   CALL timing_stop('dia_ar5')
      !
   END SUBROUTINE dia_ar5


   SUBROUTINE dia_ar5_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dia_ar5_init  ***
      !!                   
      !! ** Purpose :   initialization for AR5 diagnostic computation
      !!----------------------------------------------------------------------
      INTEGER  ::   inum
      INTEGER  ::   ik, idep
      INTEGER  ::   ji, jj, jk  ! dummy loop indices
      REAL(wp) ::   zztmp  
      REAL(wp), POINTER, DIMENSION(:,:,:,:) ::   zsaldta   ! Jan/Dec levitus salinity
      REAL(wp), POINTER, DIMENSION(:,:)     ::   zvol0    
      ! reading initial file
      LOGICAL  ::   ln_tsd_init      !: T & S data flag
      LOGICAL  ::   ln_tsd_tradmp    !: internal damping toward input data flag
      CHARACTER(len=100)            ::   cn_dir
      TYPE(FLD_N)                   ::  sn_tem,sn_sal
      INTEGER  ::   ios=0

      NAMELIST/namtsd/ ln_tsd_init,ln_tsd_tradmp,cn_dir,sn_tem,sn_sal
      !

      REWIND( numnam_ref )              ! Namelist namtsd in reference namelist :
      READ  ( numnam_ref, namtsd, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , ' namtsd in reference namelist for dia_ar5', lwp )
      REWIND( numnam_cfg )              ! Namelist namtsd in configuration namelist : Parameters of the run
      READ  ( numnam_cfg, namtsd, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , ' namtsd in configuration namelist for dia_ar5', lwp )
      IF(lwm) WRITE ( numond, namtsd )
      !
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )   CALL timing_start('dia_ar5_init')
      !
      CALL wrk_alloc( jpi, jpj, jpk, 2, zsaldta )
      CALL wrk_alloc( jpi, jpj, zvol0 )
      !                                      ! allocate dia_ar5 arrays
      IF( dia_ar5_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'dia_ar5_init : unable to allocate arrays' )

      area(:,:) = e1t(:,:) * e2t(:,:) 
      area_tot  = glob_sum( area(:,:) )

      zvol0 (:,:) = 0._wp
      thick0(:,:) = 0._wp
      DO jk = 1, jpkm1
         DO jj = 1, jpj               ! interpolation of salinity at the last ocean level (i.e. the partial step)
            DO ji = 1, jpi
               idep = tmask(ji,jj,jk) * e3t_0(ji,jj,jk)
               zvol0 (ji,jj) = zvol0 (ji,jj) +  idep * area(ji,jj)
               thick0(ji,jj) = thick0(ji,jj) +  idep    
            END DO
         END DO
      END DO
      vol0 = glob_sum( zvol0 )

      IF( iom_use('sshthster')) THEN
         CALL iom_open ( 'sali_ref_clim_monthly', inum )
         CALL iom_get  ( inum, jpdom_data, 'vosaline' , zsaldta(:,:,:,1), 1  )
         CALL iom_get  ( inum, jpdom_data, 'vosaline' , zsaldta(:,:,:,2), 12 )
         CALL iom_close( inum )

         sn0(:,:,:) = 0.5_wp * ( zsaldta(:,:,:,1) + zsaldta(:,:,:,2) )        
         sn0(:,:,:) = sn0(:,:,:) * tmask(:,:,:)
         IF( ln_zps ) THEN               ! z-coord. partial steps
            DO jj = 1, jpj               ! interpolation of salinity at the last ocean level (i.e. the partial step)
               DO ji = 1, jpi
                  ik = mbkt(ji,jj)
                  IF( ik > 1 ) THEN
                     zztmp = ( gdept_1d(ik) - gdept_0(ji,jj,ik) ) / ( gdept_1d(ik) - gdept_1d(ik-1) )
                     sn0(ji,jj,ik) = ( 1._wp - zztmp ) * sn0(ji,jj,ik) + zztmp * sn0(ji,jj,ik-1)
                  ENDIF
               END DO
            END DO
         ENDIF
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj, jpk, jpts, zsaldta )
      CALL wrk_dealloc( jpi, jpj, zvol0 )
      !
      IF( nn_timing == 1 )   CALL timing_stop('dia_ar5_init')
      !
   END SUBROUTINE dia_ar5_init

#else
   !!----------------------------------------------------------------------
   !!   Default option :                                         NO diaar5
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER :: lk_diaar5 = .FALSE.   ! coupled flag
CONTAINS
   SUBROUTINE dia_ar5_init    ! Dummy routine
   END SUBROUTINE dia_ar5_init
   SUBROUTINE dia_ar5( kt )   ! Empty routine
      INTEGER ::   kt
      WRITE(*,*) 'dia_ar5: You should not have seen this print! error?', kt
   END SUBROUTINE dia_ar5
#endif

   !!======================================================================
END MODULE diaar5
