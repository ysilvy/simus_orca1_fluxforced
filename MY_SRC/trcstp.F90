MODULE trcstp
   !!======================================================================
   !!                       ***  MODULE trcstp  ***
   !! Time-stepping    : time loop of opa for passive tracer
   !!======================================================================
   !! History :  1.0  !  2004-03  (C. Ethe)  Original
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   trc_stp      : passive tracer system time-stepping
   !!----------------------------------------------------------------------
   USE oce_trc          ! ocean dynamics and active tracers variables
   USE sbc_oce
   USE trc
   USE trctrp           ! passive tracers transport
   USE trcsms           ! passive tracers sources and sinks
   USE prtctl_trc       ! Print control for debbuging
   USE trcdia
   USE trcwri
   USE trcrst
   USE trdtrc_oce
   USE trdmxl_trc
   USE iom
   USE in_out_manager
   USE trcsub

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_stp    ! called by step

   REAL(wp), DIMENSION(:,:,:), SAVE, ALLOCATABLE ::   qsr_arr ! save qsr during TOP time-step
   REAL(wp) :: rdt_sampl
   INTEGER  :: nb_rec_per_day, ktdcy
   REAL(wp) :: rsecfst, rseclast
   LOGICAL  :: llnew

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcstp.F90 7654 2017-02-07 16:16:04Z cetlod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_stp( kt, kindic )
      !!-------------------------------------------------------------------
      !!                     ***  ROUTINE trc_stp  ***
      !!                      
      !! ** Purpose : Time loop of opa for passive tracer
      !! 
      !! ** Method  : 
      !!              Compute the passive tracers trends 
      !!              Update the passive tracers
      !!-------------------------------------------------------------------
      INTEGER, INTENT( in ) ::  kt      ! ocean time-step index
      INTEGER, INTENT( inout ) ::   kindic  ! indicator 
      INTEGER               ::  jk, jn  ! dummy loop indices
      REAL(wp)              ::  ztrai
      CHARACTER (len=25)    ::  charout 
      !!-------------------------------------------------------------------
      !
      IF( nn_timing == 1 )   CALL timing_start('trc_stp')
      !
      IF( kt == nittrc000 .AND. lk_trdmxl_trc )  CALL trd_mxl_trc_init    ! trends: Mixed-layer
      !
      IF( lk_vvl ) THEN                                                   ! update ocean volume due to ssh temporal evolution
         DO jk = 1, jpk
            cvol(:,:,jk) = e1e2t(:,:) * fse3t(:,:,jk) * tmask(:,:,jk)
         END DO
         IF( lk_degrad )  cvol(:,:,:) = cvol(:,:,:) * facvol(:,:,:)       ! degrad option: reduction by facvol
         areatot         = glob_sum( cvol(:,:,:) )
      ENDIF
      !
      IF( l_trcdm2dc )   CALL trc_mean_qsr( kt )
      !    
      IF( nn_dttrc /= 1 )   CALL trc_sub_stp( kt )  ! averaging physical variables for sub-stepping
      !    
      IF( MOD( kt , nn_dttrc ) == 0 ) THEN      ! only every nn_dttrc time step
         !
         IF(ln_ctl) THEN
            WRITE(charout,FMT="('kt =', I4,'  d/m/y =',I2,I2,I4)") kt, nday, nmonth, nyear
            CALL prt_ctl_trc_info(charout)
         ENDIF
         !
         tra(:,:,:,:) = 0.e0
         !
         IF( .NOT.lk_offline )     CALL trc_rst_opn  ( kt )       ! Open tracer restart file 
         IF( lrst_trc )            CALL trc_rst_cal  ( kt, 'WRITE' )   ! calendar
         IF( lk_iomput ) THEN  ;   CALL trc_wri      ( kt )       ! output of passive tracers with iom I/O manager
         ELSE                  ;   CALL trc_dia      ( kt )       ! output of passive tracers with old I/O manager
         ENDIF
                                   CALL trc_sms      ( kt )       ! tracers: sinks and sources
                                   CALL trc_trp      ( kt )       ! transport of passive tracers

         IF( kt == nittrc000 ) THEN
            CALL iom_close( numrtr )       ! close input tracer restart file
            IF(lwm) CALL FLUSH( numont )   ! flush namelist output
         ENDIF
         IF( lrst_trc )            CALL trc_rst_wri  ( kt )       ! write tracer restart file
         IF( lk_trdmxl_trc  )      CALL trd_mxl_trc  ( kt )       ! trends: Mixed-layer
         !
         IF( nn_dttrc /= 1   )     CALL trc_sub_reset( kt )       ! resetting physical variables when sub-stepping
         !
      ENDIF
      !
      CALL trc_stp_ctl( kt, kindic )
      IF( lk_iomput .AND. kindic < 0 )  THEN
         CALL trc_wri_state( 'output.trc.abort', kt )
         nstop = nstop + 1
      ENDIF

      !
      ztrai = 0._wp                                                   !  content of all tracers
      DO jn = 1, jptra
         ztrai = ztrai + glob_sum( trn(:,:,:,jn) * cvol(:,:,:)   )
      END DO
      IF( lwp ) WRITE(numstr,9300) kt,  ztrai / areatot
9300  FORMAT(i10,D23.16)
      !
      IF( nn_timing == 1 )   CALL timing_stop('trc_stp')
      !
   END SUBROUTINE trc_stp

   SUBROUTINE trc_mean_qsr( kt )
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE trc_mean_qsr  ***
      !!
      !! ** Purpose :  Compute daily mean qsr for biogeochemical model in case
      !!               of diurnal cycle
      !!
      !! ** Method  : store in TOP the qsr every hour ( or every time-step if the latter 
      !!              is greater than 1 hour ) and then, compute the  mean with 
      !!              a moving average over 24 hours. 
      !!              In coupled mode, the sampling is done at every coupling frequency 
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt
      INTEGER  :: jn
      REAL(wp) :: zkt, zrec
      CHARACTER(len=1)               ::   cl1                      ! 1 character
      CHARACTER(len=2)               ::   cl2                      ! 2 characters

      IF( kt == nittrc000 ) THEN
         IF( ln_cpl )  THEN  
            rdt_sampl = rday / ncpl_qsr_freq
            nb_rec_per_day = ncpl_qsr_freq
         ELSE  
            rdt_sampl = MAX( 3600., rdttrc(1) )
            nb_rec_per_day = INT( rday / rdt_sampl )
         ENDIF
         !
         IF( lwp ) THEN
            WRITE(numout,*) 
            WRITE(numout,*) ' Sampling frequency dt = ', rdt_sampl, 's','   Number of sampling per day  nrec = ', nb_rec_per_day
            WRITE(numout,*) 
         ENDIF
         !
         ALLOCATE( qsr_arr(jpi,jpj,nb_rec_per_day ) )
         !
         !                                            !* Restart: read in restart file
         IF( ln_rsttr .AND. nn_rsttr /= 0 .AND. iom_varid( numrtr, 'qsr_mean' , ldstop = .FALSE. ) > 0  &
           &                              .AND. iom_varid( numrtr, 'qsr_arr_1', ldstop = .FALSE. ) > 0  &
           &                              .AND. iom_varid( numrtr, 'ktdcy'    , ldstop = .FALSE. ) > 0  &
           &                              .AND. iom_varid( numrtr, 'nrdcy'    , ldstop = .FALSE. ) > 0  ) THEN

            CALL iom_get( numrtr, 'ktdcy', zkt )  
            rsecfst = INT( zkt ) * rdttrc(1)
            IF(lwp) WRITE(numout,*) 'trc_qsr_mean:   qsr_mean read in the restart file at time-step rsecfst =', rsecfst, ' s '
            CALL iom_get( numrtr, jpdom_autoglo, 'qsr_mean', qsr_mean )   !  A mean of qsr
            CALL iom_get( numrtr, 'nrdcy', zrec )   !  Number of record per days
            IF( INT( zrec ) == nb_rec_per_day ) THEN
               DO jn = 1, nb_rec_per_day 
                  IF( jn <= 9 )  THEN
                    WRITE(cl1,'(i1)') jn
                    CALL iom_get( numrtr, jpdom_autoglo, 'qsr_arr_'//cl1, qsr_arr(:,:,jn) )   !  A mean of qsr
                  ELSE
                    WRITE(cl2,'(i2.2)') jn
                    CALL iom_get( numrtr, jpdom_autoglo, 'qsr_arr_'//cl2, qsr_arr(:,:,jn) )   !  A mean of qsr
                  ENDIF
              ENDDO
            ELSE
               DO jn = 1, nb_rec_per_day
                  qsr_arr(:,:,jn) = qsr_mean(:,:)
               ENDDO
            ENDIF
         ELSE                                         !* no restart: set from nit000 values
            IF(lwp) WRITE(numout,*) 'trc_qsr_mean:   qsr_mean set to nit000 values'
            rsecfst  = kt * rdttrc(1)
            !
            qsr_mean(:,:) = qsr(:,:)
            DO jn = 1, nb_rec_per_day
               qsr_arr(:,:,jn) = qsr_mean(:,:)
            ENDDO
         ENDIF
         !
      ENDIF
      !
      rseclast = kt * rdttrc(1)
      !
      llnew   = ( rseclast - rsecfst ) .ge.  rdt_sampl    !   new shortwave to store
      IF( llnew ) THEN
          ktdcy = kt
          IF( lwp .AND. kt < nittrc000 + 100 ) WRITE(numout,*) ' New shortwave to sample for TOP at time kt = ', ktdcy, &
             &                      ' time = ', rseclast/3600.,'hours '
          rsecfst = rseclast
          DO jn = 1, nb_rec_per_day - 1
             qsr_arr(:,:,jn) = qsr_arr(:,:,jn+1)
          ENDDO
          qsr_arr (:,:,nb_rec_per_day) = qsr(:,:)
          qsr_mean(:,:                ) = SUM( qsr_arr(:,:,:), 3 ) / nb_rec_per_day
      ENDIF
      !
      IF( lrst_trc ) THEN    !* Write the mean of qsr in restart file 
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'trc_mean_qsr : write qsr_mean in restart file  kt =', kt
         IF(lwp) WRITE(numout,*) '~~~~~~~'
         zkt  = REAL( ktdcy, wp )
         zrec = REAL( nb_rec_per_day, wp )
         CALL iom_rstput( kt, nitrst, numrtw, 'ktdcy', zkt  )
         CALL iom_rstput( kt, nitrst, numrtw, 'nrdcy', zrec )
          DO jn = 1, nb_rec_per_day 
             IF( jn <= 9 )  THEN
               WRITE(cl1,'(i1)') jn
               CALL iom_rstput( kt, nitrst, numrtw, 'qsr_arr_'//cl1, qsr_arr(:,:,jn) )
             ELSE
               WRITE(cl2,'(i2.2)') jn
               CALL iom_rstput( kt, nitrst, numrtw, 'qsr_arr_'//cl2, qsr_arr(:,:,jn) )
             ENDIF
         ENDDO
         CALL iom_rstput( kt, nitrst, numrtw, 'qsr_mean', qsr_mean(:,:) )
      ENDIF
      !
   END SUBROUTINE trc_mean_qsr

   SUBROUTINE trc_stp_ctl( kt, kindic )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE stp_ctl  ***
      !!                     
      !! ** Purpose :   Control the run
      !!
      !! ** Method  : - Save the time step in numstp
      !!              - Print it each 50 time steps
      !!              - Print solver statistics in numsol 
      !!              - Stop the run IF problem for the solver ( indec < 0 )
      !!
      !! ** Actions :   'time.step' file containing the last ocean time-step
      !!                
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step index
      INTEGER, INTENT( inout ) ::   kindic  ! indicator of solver convergence
      !!
      INTEGER  ::   ji, jj, jk, jn              ! dummy loop indices
      INTEGER  ::   ii, ij, ik, itrc           ! temporary integers
      INTEGER, DIMENSION(3) ::   iloc      ! 
      REAL(wp) ::   zmax, zhuge     ! temporary scalars
      REAL(wp), DIMENSION(jptra) ::   ztrc     ! temporary scalars
      !!----------------------------------------------------------------------


      DO jn = 1, jptra
         ztrc(jn)  = MAXVAL( trn(:,:,:,jn), mask= ((tmask*SPREAD(tmask_i,DIM=3,NCOPIES=jpk).NE.0.)) )
         IF( lk_mpp ) CALL mpp_max( ztrc(jn) )      ! max over the global domain
      END DO

      zmax = ztrc(1)
      itrc = 1
      DO jn = 2, jptra
         IF( ztrc(jn) > zmax ) THEN
            zmax = ztrc(jn) 
            itrc = jn
         ENDIF
      END DO 

      zhuge = 1.e+10
      IF( ( zmax > zhuge ) .OR. isnan( zmax )  ) THEN
         IF( lk_mpp ) THEN
            CALL mpp_maxloc( trn(:,:,:,itrc),tmask, zmax,ii,ij,ik)
         ELSE
            iloc = MAXLOC( trn(:,:,:,itrc) )
            ii = iloc(1) + nimpp - 1
            ij = iloc(2) + njmpp - 1
            ik = iloc(3)
         ENDIF
         IF(lwp) THEN
            WRITE(numout,cform_err)
            WRITE(numout,*) "  trc_stp_ctl : passive tracer  ", TRIM( ctrcnm(itrc) ), "  is too big.  "
            WRITE(numout,*) ' ============ '
            WRITE(numout,9400) kt, zmax, ii, ij, ik
            WRITE(numout,*)
            WRITE(numout,*) '          output of last fields in numwso'
         ENDIF
         kindic = -3
      ENDIF
9400  FORMAT (' kt=',i6,' trcmax= ',D20.8,', i j k: ',3i5)
   END SUBROUTINE trc_stp_ctl

#else
   !!----------------------------------------------------------------------
   !!   Default key                                     NO passive tracers
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_stp( kt )        ! Empty routine
      WRITE(*,*) 'trc_stp: You should not have seen this print! error?', kt
   END SUBROUTINE trc_stp
#endif

   !!======================================================================
END MODULE trcstp
