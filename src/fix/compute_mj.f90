FUNCTION compute_mj(j,l,m)
   !-----------------------------------------------------------------------
   USE kinds, ONLY: DP
   IMPLICIT NONE
   !
   REAL(DP) :: compute_mj, j
   INTEGER  :: l, m

   IF (abs(j-l-0.5d0)<1.d-4) THEN
       compute_mj=m+0.5d0
   ELSEIF (abs(j-l+0.5d0)<1.d-4) THEN
      compute_mj=m-0.5d0
   ELSE
      CALL errore('compute_mj','l and j not compatible',1)
   ENDIF

   RETURN
END FUNCTION compute_mj
