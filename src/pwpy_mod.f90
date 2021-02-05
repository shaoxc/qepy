MODULE pwpy_mod
   USE kinds,                ONLY : DP
   USE pwpy_scatter_mod, ONLY : gather_grid, scatter_grid
   IMPLICIT NONE
   PUBLIC

   type, public :: embed_base
      real(kind=dp), allocatable      :: extpot(:)
      real(kind=dp)                   :: extene = 0.0
      integer                         :: exttype = 0
      logical                         :: initial = .true.
      real(kind=dp)                   :: mix_coef = -1.0
      logical                         :: finish = .false.
      real(kind=dp)                   :: etotal = 0.0
      real(kind=dp)                   :: dnorm = 1.0
      logical                         :: lewald = .true.
      logical                         :: nlpp = .true.
      real(kind=dp)                   :: diag_conv = 1.D-2
   end type embed_base

CONTAINS
   SUBROUTINE pwpy_init_pointer()
      use scf, only: rho,v,vnew
      !
      IMPLICIT NONE
      REAL(DP), POINTER :: scf__rho__of_r(:,:)
      REAL(DP), POINTER :: scf__v__of_r(:,:)
      REAL(DP), POINTER :: scf__vnew__of_r(:,:)
      !
      IF (ALLOCATED(rho%of_r)) THEN
         associate(scf__rho__of_r => rho%of_r)
         end associate
      ENDIF
      IF (ALLOCATED(v%of_r)) THEN
         associate(scf__v__of_r => v%of_r)
         end associate
      ENDIF
      IF (ALLOCATED(vnew%of_r)) THEN
         associate(scf__vnew__of_r => vnew%of_r)
         end associate
      ENDIF
   END SUBROUTINE

   SUBROUTINE pwpy_get_rho(rhor)
      USE kinds,                ONLY : DP
      use scf, only: rho !! the charge density and its other components
      USE fft_base,         ONLY : dfftp, dffts
      USE lsda_mod,           ONLY : lsda, nspin, current_spin, isk
      !
      IMPLICIT NONE
      REAL(DP), INTENT(OUT) :: rhor(:,:)
      !REAL(DP), INTENT(OUT) :: rhor(dfftp%nr1x * dfftp%nr2x * dfftp%nr3x, nspin)
      !
      INTEGER  :: ispin
      !
      !print *, 'get_rho_IN',minval(rho%of_r),maxval(rho%of_r),sum(rho%of_r)
      DO ispin = 1, nspin
         CALL gather_grid(dfftp, rho%of_r(:,ispin), rhor(:,ispin))
      END DO
      !print *, 'get_rho_OUT',minval(rhor),maxval(rhor),sum(rhor)
   END SUBROUTINE

   SUBROUTINE pwpy_set_rho(rhor)
      USE kinds,                ONLY : DP
      USE fft_rho,              ONLY : rho_g2r, rho_r2g
      USE fft_base,         ONLY : dfftp, dffts
      use scf, only: rho !! the charge density and its other components
      USE lsda_mod,           ONLY : lsda, nspin, current_spin, isk
      !
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: rhor(:,:)
      !
      INTEGER  :: ispin
      DO ispin = 1, nspin
         CALL scatter_grid(dfftp, rhor(:,ispin), rho%of_r(:,ispin))
      END DO
      CALL rho_r2g(dfftp, rho%of_r, rho%of_g )
   END SUBROUTINE

   SUBROUTINE pwpy_get_rho_core(rhoc)
      USE kinds,                ONLY : DP
      use scf, only: rho_core !! the core charge in real space
      USE fft_base,         ONLY : dfftp, dffts
      IMPLICIT NONE
      REAL(DP), INTENT(OUT) :: rhoc(:)
      !
      CALL gather_grid(dfftp, rho_core, rhoc)
   END SUBROUTINE

   SUBROUTINE pwpy_set_rho_core(rhoc)
      USE kinds,                ONLY : DP
      use scf, only: rho_core !! the core charge in real space
      USE fft_base,         ONLY : dfftp, dffts
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: rhoc(:)
      !
      CALL scatter_grid(dfftp, rhoc, rho_core)
   END SUBROUTINE

   SUBROUTINE pwpy_set_extpot(embed, vin)
      USE kinds,                ONLY : DP
      USE fft_rho,              ONLY : rho_g2r, rho_r2g
      USE fft_base,         ONLY : dfftp, dffts
      !
      IMPLICIT NONE
      TYPE(embed_base), INTENT(INOUT) :: embed
      REAL(DP), INTENT(IN) :: vin(:)
      !
      IF (ALLOCATED(embed%extpot)) THEN
         IF (SIZE(embed%extpot) /= dfftp%nnr) DEALLOCATE(embed%extpot)
      ENDIF
      IF (.NOT.ALLOCATED(embed%extpot)) ALLOCATE(embed%extpot(dfftp%nnr))
      CALL scatter_grid(dfftp, vin, embed%extpot)
   END SUBROUTINE

   SUBROUTINE pwpy_get_grid(nr)
      USE kinds,                ONLY : DP
      USE fft_base,         ONLY : dfftp, dffts
      !
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: nr(3)
      !
      nr=(/dfftp%nr1, dfftp%nr2, dfftp%nr3/)
   END SUBROUTINE

   SUBROUTINE pwpy_set_stdout(fname, uni)
      USE io_global,     ONLY : stdout, ionode
      !
      INTEGER                  :: ierr
      CHARACTER(LEN=*),INTENT(IN),OPTIONAL  :: fname
      INTEGER,INTENT(in),OPTIONAL :: uni
      !
      IF ( .NOT. present(fname) ) return
      IF ( present(uni) ) THEN
         stdout=uni
      ELSE
         stdout=666
      ENDIF

      IF(ionode) THEN
         OPEN (UNIT = stdout, FILE = TRIM(fname), FORM = 'formatted', STATUS = 'unknown', iostat = ierr )
      ENDIF
   END SUBROUTINE

   SUBROUTINE pwpy_write_stdout(fstr)
      USE io_global,     ONLY : stdout, ionode
      !
      INTEGER                  :: ierr
      CHARACTER(LEN=*),INTENT(IN)  :: fstr
      !
      IF(ionode) WRITE(stdout,'(A)') fstr
   END SUBROUTINE

   SUBROUTINE pwpy_close_stdout(fname)
      USE io_global,     ONLY : stdout, ionode
      !
      INTEGER                  :: ierr
      CHARACTER(LEN=*),INTENT(IN)  :: fname
      !
      IF(ionode) close(stdout)
   END SUBROUTINE

   SUBROUTINE pwpy_update_ions(embed, pos)
      ! This is function Combined 'run_pwscf' and 'move_ions'
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
      !
      CALL update_file()
      IF ( ionode ) THEN
         tau(:,:)=pos(:,:)
         CALL checkallsym( nat, tau, ityp)
      ENDIF
      CALL mp_bcast( tau, ionode_id, intra_image_comm )
      CALL punch( 'config-nowf' )
      !
      IF ( treinit_gvecs ) THEN
         CALL reset_gvectors()
      ELSE
         CALL update_pot()
         CALL pwpy_hinit1(embed%exttype)
      END IF
   END SUBROUTINE

   SUBROUTINE pwpy_set_mod_float(mod_name, param, value)
      !
      CHARACTER(LEN=*),INTENT(IN)  :: mod_name
      CHARACTER(LEN=*),INTENT(IN)  :: param
      REAL(DP),INTENT(IN)  :: value
      !
   END SUBROUTINE

END MODULE pwpy_mod
