MODULE trcnam_age
   !!======================================================================
   !!                         ***  MODULE trcnam_age  ***
   !! TOP :   initialisation of some run parameters for Age tracer
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) 
   !!----------------------------------------------------------------------
#if defined key_age
   !!----------------------------------------------------------------------
   !!   'key_age'                                               AGE tracers
   !!----------------------------------------------------------------------
   !! trc_nam_age      : AGE  tracer initialisation
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE trcsms_age      ! AGE specific variable
   USE trc

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_age   ! called by trcnam.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcnam_age.F90 8353 2017-07-19 14:41:00Z lovato $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nam_age
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE trc_nam_age  ***
      !!                 
      !! ** Purpose :   Definition some run parameter for AGE model
      !!
      !! ** input   :   Namelist namage
      !!----------------------------------------------------------------------
      ! Variable setting
      ctrcnm    (jp_age0) = 'Age'
      ctrcln    (jp_age0) = 'Sea water age since surface contact'
      ctrcun    (jp_age0) = 'year'
      ln_trc_ini(jp_age0) = .false.

      rn_age_depth      = 10
      rn_age_kill_rate  = -1./7200

      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) ' trc_nam_age: Read namage, namelist for Age passive tracer'
         WRITE(numout,*) ' ~~~~~~~'
         WRITE(numout,*) '  depth over which age tracer reset to zero                              rn_age_depth      = ', rn_age_depth 
         WRITE(numout,*) '  recip of relax. timescale (s) for age tracer shallower than age_depth  rn_age_kill_rate  = ', rn_age_kill_rate 
      ENDIF

   END SUBROUTINE trc_nam_age
   
#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                                No AGE
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_nam_age                      ! Empty routine
   END  SUBROUTINE  trc_nam_age
#endif  

   !!======================================================================
END MODULE trcnam_age
