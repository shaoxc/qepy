MODULE qepy_tddft_mod
   USE kinds,                ONLY : DP
   !
   IMPLICIT NONE
   PUBLIC
   !
CONTAINS

   SUBROUTINE qepy_cetddft_wfc2rho(iunit)
      USE kinds,                ONLY : DP
      USE io_files,             ONLY : iunwfc, nwordwfc
      USE wavefunctions,        ONLY : evc
      USE klist,                ONLY : nks
      USE buffers,              ONLY : get_buffer
      !
      USE tddft_module, ONLY: iunevcn, l_tddft_restart
      IMPLICIT NONE
      INTEGER, INTENT(IN),OPTIONAL :: iunit
      !
      IF ( present(iunit) ) THEN
        iunwfc = iunit
      ELSE
        iunwfc = iunevcn
     ENDIF
     if (nks==1) then
        call get_buffer(evc, nwordwfc, iunwfc, 1)
     endif
     call sum_band()
     l_tddft_restart = .true.
  END SUBROUTINE

END MODULE qepy_tddft_mod
