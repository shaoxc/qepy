MODULE qepy_api
   USE kinds,                ONLY : DP
   USE qepy_common, ONLY : embed_base
   !
   IMPLICIT NONE
   PUBLIC
   !
CONTAINS

   SUBROUTINE qepy_update_ions(embed, pos, ikind)
      !-----------------------------------------------------------------------
      ! This is function Combined 'run_pwscf' and 'move_ions'
      ! ikind = 0  all
      ! ikind = 1  atomic configuration dependent information
      !-----------------------------------------------------------------------
      USE mp_images,            ONLY : intra_image_comm
      USE extrapolation,        ONLY : update_file, update_pot
      USE io_global,            ONLY : ionode_id, ionode
      USE ions_base,            ONLY : nat, ityp, tau
      USE symm_base,            ONLY : checkallsym
      USE mp,                   ONLY : mp_bcast
      USE control_flags,        ONLY : treinit_gvecs
      !
      INTEGER                  :: ierr
      TYPE(embed_base), INTENT(INOUT) :: embed
      REAL(DP), INTENT(IN) :: pos(:,:)
      INTEGER,INTENT(IN),OPTIONAL  :: ikind
      !
      INTEGER   :: iflag
      !
      IF ( present(ikind) ) THEN
         iflag = ikind
      ELSE
         iflag = 0
      ENDIF

      CALL update_file()

      IF ( ionode ) THEN
         tau(:,:)=pos(:,:)
         CALL checkallsym( nat, tau, ityp)
      ENDIF
      CALL mp_bcast( tau, ionode_id, intra_image_comm )
      IF (iflag == 0 ) THEN
         CALL punch( 'config-nowf' )
         IF ( treinit_gvecs ) THEN
            CALL reset_gvectors()
         ELSE
            CALL update_pot()
            CALL qepy_hinit1(embed%exttype)
         END IF
      ELSE IF (iflag == 1 ) THEN
         CALL set_rhoc()
         CALL qepy_hinit1(embed%exttype)
      END IF
   END SUBROUTINE

END MODULE qepy_api
