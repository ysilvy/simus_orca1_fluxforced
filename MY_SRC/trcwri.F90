MODULE trcwri
   !!======================================================================
   !!                       *** MODULE trcwri ***
   !!    TOP :   Output of passive tracers
   !!======================================================================
   !! History :   1.0  !  2009-05 (C. Ethe)  Original code
   !!----------------------------------------------------------------------
#if defined key_top && defined key_iomput
   !!----------------------------------------------------------------------
   !!   'key_top'                                           TOP models
   !!----------------------------------------------------------------------
   !! trc_wri_trc   :  outputs of concentration fields
   !!----------------------------------------------------------------------
   USE dom_oce     ! ocean space and time domain variables
   USE oce_trc     ! shared variables between ocean and passive tracers
   USE trc         ! passive tracers common variables 
   USE iom         ! I/O manager
   USE ioipsl
   USE dianam      ! Output file name
   USE trcwri_pisces
   USE trcwri_cfc
   USE trcwri_c14b
   USE trcwri_age
   USE trcwri_my_trc

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_wri      
   PUBLIC trc_wri_state      

   !! * Substitutions
#  include "top_substitute.h90"

CONTAINS

   SUBROUTINE trc_wri( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_wri  ***
      !! 
      !! ** Purpose :   output passive tracers fields and dynamical trends
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in )     :: kt
      !
      INTEGER                   :: jn
      CHARACTER (len=20)        :: cltra
      CHARACTER (len=40)        :: clhstnam
      INTEGER ::   inum = 11            ! temporary logical unit
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_wri')
      !
      IF( lk_offline .AND. kt == nittrc000 .AND. lwp ) THEN    ! WRITE root name in date.file for use by postpro
         CALL dia_nam( clhstnam, nn_writetrc,' ' )
         CALL ctl_opn( inum, 'date.file', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, narea )
         WRITE(inum,*) clhstnam
         CLOSE(inum)
      ENDIF
      ! write the tracer concentrations in the file
      ! ---------------------------------------
      IF( lk_pisces  )   CALL trc_wri_pisces     ! PISCES 
      IF( lk_my_trc  )   CALL trc_wri_my_trc     ! MY_TRC  tracers
      IF( lk_cfc     )   CALL trc_wri_cfc        ! surface fluxes of CFC
      IF( lk_c14b    )   CALL trc_wri_c14b       ! surface fluxes of C14
      IF( lk_age     )   CALL trc_wri_age        ! AGE tracer
      !
      IF( nn_timing == 1 )  CALL timing_stop('trc_wri')
      !
   END SUBROUTINE trc_wri

   SUBROUTINE trc_wri_state( cdfile_name, kt )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE dia_wri_state  ***
      !!        
      !! ** Purpose :   create a NetCDF file named cdfile_name which contains 
      !!      the instantaneous ocean state and forcing fields.
      !!        Used to find errors in the initial state or save the last
      !!      ocean state in case of abnormal end of a simulation
      !!
      !! ** Method  :   NetCDF files using ioipsl
      !!      File 'output.init.nc'  is created if ninist = 1 (namelist)
      !!      File 'output.abort.trc.nc' is created in case of abnormal job end
      !!----------------------------------------------------------------------
      CHARACTER (len=* ), INTENT( in ) ::   cdfile_name      ! name of the file created
      INTEGER           , INTENT( in ) ::   kt               ! ocean time-step index
      !! 
      CHARACTER (len=32) :: clname
      CHARACTER (len=40) :: clop
      CHARACTER (len=20) :: cltra, cltrau
      CHARACTER (len=80) :: cltral
      INTEGER  ::   id_i , nz_i, nh_i, jn       
      INTEGER, DIMENSION(1) ::   idex             ! local workspace
      REAL(wp) ::   zsto, zout, zmax, zjulian, zdt
      !!----------------------------------------------------------------------
      ! 
      ! -----------------

      ! Define name, frequency of output and means
      clname = cdfile_name
      zdt  = rdt 
      zsto = rdt 
      clop = "inst(x)"           ! no use of the mask value (require less cpu time)
      zout = rdt

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'trc_wri_state : single instantaneous ocean state'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~   and forcing fields file created '
      IF(lwp) WRITE(numout,*) '                and named :', clname, '.nc'

      ! Compute julian date from starting date of the run
      CALL ymds2ju( nyear, nmonth, nday, rdt, zjulian )
      zjulian = zjulian - adatrj   !   set calendar origin to the beginning of the experiment

      CALL histbeg( clname, jpi, glamt, jpj, gphit,   &
          1, jpi, 1, jpj, nit000-1, zjulian, zdt, nh_i, id_i, domain_id=nidom, snc4chunks=snc4set ) ! Horizontal grid : glamt and gphit
      CALL histvert( id_i, "deptht", "Vertical T levels", "m", jpk, gdept_1d, nz_i, "down")  ! Vertical grid : gdept



      ! Declare all the output fields as NETCDF variables
      DO jn = 1, jptra
         cltra  = TRIM( ctrcnm(jn) )   ! short title for tracer
         cltral = TRIM( ctrcln(jn) )   ! long title for tracer
         cltrau = TRIM( ctrcun(jn) )   ! UNIT for tracer
         CALL histdef( id_i, cltra, cltral, cltrau, jpi, jpj, nh_i,  &
            &          jpk, 1, jpk,  nz_i, 32, clop, zsto, zout ) 
      END DO

      CALL histend( id_i, snc4chunks=snc4set )


      DO jn = 1, jptra
         cltra  = TRIM( ctrcnm(jn) )   ! short title for tracer
         CALL histwrite( id_i, cltra, kt, trn(:,:,:,jn), jpi*jpj*jpk, idex )
      END DO

      CALL histclo( id_i )
   END SUBROUTINE trc_wri_state

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                     No passive tracer
   !!----------------------------------------------------------------------
   PUBLIC trc_wri
CONTAINS
   SUBROUTINE trc_wri( kt )                     ! Empty routine   
   INTEGER, INTENT(in) :: kt
   END SUBROUTINE trc_wri

   SUBROUTINE trc_wri_state( cd_name, kt )                     ! Empty routine   
   CHARACTER (len=* ), INTENT( in ) ::   cd_name      ! name of the file created
   INTEGER, INTENT(in) :: kt
   END SUBROUTINE trc_wri_state
#endif

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcwri.F90 8353 2017-07-19 14:41:00Z lovato $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE trcwri
