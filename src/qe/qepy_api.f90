MODULE qepy_api
   USE kinds,                ONLY : DP
   USE qepy_common, ONLY : embed_base
   !
   IMPLICIT NONE
   PUBLIC
   !
CONTAINS

   SUBROUTINE qepy_update_ions(embed, pos, ikind, lattice)
      !-----------------------------------------------------------------------
      ! This is function Combined 'run_pwscf' and 'move_ions'.
      !***********************************************************************
      ! pos:
      !   ionic positions in bohr
      ! ikind:
      !   ikind = 0  all
      !   ikind = 1  atomic configuration dependent information
      ! lattice:
      !   lattice parameter in bohr
      !***********************************************************************
      !-----------------------------------------------------------------------
      USE mp_images,            ONLY : intra_image_comm
      USE extrapolation,        ONLY : update_file, update_pot
      USE io_global,            ONLY : ionode_id, ionode
      USE ions_base,            ONLY : nat, ityp, tau
      USE symm_base,            ONLY : checkallsym
      USE mp,                   ONLY : mp_bcast
      USE control_flags,        ONLY : treinit_gvecs
      USE cell_base,            ONLY : alat, at, bg, omega, cell_force, &
                                     fix_volume, fix_area, ibrav, enforce_ibrav
      USE cellmd,               ONLY : omega_old, at_old, press, lmovecell, calc, cell_factor
      !
      INTEGER                  :: ierr
      TYPE(embed_base), INTENT(INOUT) :: embed
      REAL(DP), INTENT(IN) :: pos(:,:)
      INTEGER,INTENT(IN),OPTIONAL  :: ikind
      REAL(DP), INTENT(IN), OPTIONAL  :: lattice(3,3)
      !
      INTEGER   :: iflag
      LOGICAL :: lmovecell_
      !
      IF ( present(ikind) ) THEN
         iflag = ikind
      ELSE
         iflag = 0
      ENDIF
      !
      IF ( present(lattice) ) THEN
         IF (.not. lmovecell) THEN
            call errore("qepy_update_ions","lattice update only works for calculation= 'vc-relax' and 'vc-md'.",1)
            ! This is due to the `gshells` function in the initialization will set `ngl`, which is one of the shape of `vloc`.
         ENDIF
         lmovecell_ = .TRUE.
      ELSE
         lmovecell_ = .FALSE.
      ENDIF

      CALL update_file()

      IF ( ionode ) THEN
         tau(:,:)=pos(:,:) / alat
         IF ( lmovecell_ ) THEN
            !
            IF (ALLOCATED(embed%extpot)) DEALLOCATE(embed%extpot)
            !
            at_old = at
            omega_old = omega
            IF (fix_volume) CALL impose_deviatoric_strain( alat*at, lattice )
            IF (fix_area)   CALL impose_deviatoric_strain_2d( alat*at, lattice )
            at = lattice / alat
            IF(enforce_ibrav) CALL remake_cell( ibrav, alat, at(1,1),at(1,2),at(1,3) )
            CALL recips( at(1,1),at(1,2),at(1,3), bg(1,1),bg(1,2),bg(1,3) )
            CALL volume( alat, at(1,1),at(1,2),at(1,3), omega )
            !
         ENDIF
         CALL checkallsym( nat, tau, ityp)
      ENDIF
      !
      CALL mp_bcast( tau, ionode_id, intra_image_comm )
      !
      IF ( lmovecell_ ) THEN
         !
         CALL mp_bcast( at,        ionode_id, intra_image_comm )
         CALL mp_bcast( at_old,    ionode_id, intra_image_comm )
         CALL mp_bcast( omega,     ionode_id, intra_image_comm )
         CALL mp_bcast( omega_old, ionode_id, intra_image_comm )
         CALL mp_bcast( bg,        ionode_id, intra_image_comm )
         !
      ENDIF
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
