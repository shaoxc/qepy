MODULE qepy_common
   USE kinds,                ONLY : DP
   IMPLICIT NONE
   PRIVATE
   !
   PUBLIC :: arr2pointer, allocate_extpot, allocate_extforces
   !
   type, public :: input_base
      INTEGER            :: my_world_comm = 0
      LOGICAL            :: start_images = .false.
      CHARACTER(len=256) :: filename = ''
      CHARACTER(len=256) :: code = 'QEPY'
      CHARACTER(len=256) :: tmp_dir = './'
      CHARACTER(len=256) :: wfc_dir = 'undefined'
      CHARACTER(len=256) :: prefix  = 'os'
      ! ...  If needwf=.t. performs wavefunction-related initialization as well
      LOGICAL            :: needwf
   end type input_base
   !
   type, public :: tddft_base
      logical                         :: initial = .true.
      logical                         :: finish = .false.
      integer                         :: istep = 0
      integer                         :: nstep = 1
      logical                         :: iterative = .false.
      real(kind=dp), allocatable      :: dipole(:,:)
   end type tddft_base
   !
   type, public :: embed_base
      type(input_base)                :: input
      type(tddft_base)                :: tddft
      real(kind=dp), allocatable      :: extpot(:,:)
      real(kind=dp)                   :: extene = 0.0
      integer                         :: exttype = 0
      real(kind=dp), allocatable      :: extforces(:,:)
      real(kind=dp)                   :: extstress(3,3)
      logical                         :: initial = .true.
      real(kind=dp)                   :: mix_coef = -1.0
      logical                         :: finish = .false.
      real(kind=dp)                   :: etotal = 0.0
      real(kind=dp)                   :: dnorm = 1.0
      logical                         :: lewald = .true.
      logical                         :: nlpp = .true.
      real(kind=dp)                   :: diag_conv = 1.D-2
      logical                         :: ldescf = .false.
      !! add scf correction energy
      logical                         :: iterative = .false.
      !! add correction for variational energy
      logical                         :: lmovecell = .false.
      !! allow change the cell
      logical                         :: oldxml = .false.
      !! Olderversion QE (XML file name is 'data-file.xml')
   CONTAINS
      !--------------------------------------------------------------------------------
      PROCEDURE :: allocate_extpot => allocate_extpot_class
      PROCEDURE :: allocate_extforces => allocate_extforces_class
   end type embed_base
   !
   !
   TYPE ( embed_base ), public :: messenger
   !
   !
   INTERFACE arr2pointer
      MODULE PROCEDURE arr2pointer_real_1, arr2pointer_real_2, arr2pointer_real_3, arr2pointer_real_4
   END INTERFACE
   !
CONTAINS
   !
   SUBROUTINE allocate_extpot_class(embed)
      USE kinds,                ONLY : DP
      USE fft_base,             ONLY : dfftp
      USE lsda_mod,             ONLY : lsda, nspin
      !
      IMPLICIT NONE
      CLASS(embed_base), INTENT(INOUT) :: embed
      !
      CALL allocate_extpot(embed)
      !
   END SUBROUTINE
   !
   SUBROUTINE allocate_extpot(embed)
      USE kinds,                ONLY : DP
      USE fft_base,             ONLY : dfftp
      USE lsda_mod,             ONLY : lsda, nspin
      !
      IMPLICIT NONE
      TYPE(embed_base), INTENT(INOUT) :: embed
      !
      IF (ALLOCATED(embed%extpot)) THEN
         IF (SIZE(embed%extpot, 1) /= dfftp%nnr) DEALLOCATE(embed%extpot)
      ENDIF
      IF (.NOT.ALLOCATED(embed%extpot)) THEN
         ALLOCATE(embed%extpot(dfftp%nnr, nspin))
         embed%extpot = 0.0_DP
      ENDIF
   END SUBROUTINE
   !
   SUBROUTINE allocate_extforces_class(embed)
      USE kinds,                ONLY : DP
      USE fft_base,             ONLY : dfftp
      USE ions_base,            ONLY : nat
      !
      IMPLICIT NONE
      CLASS(embed_base), INTENT(INOUT) :: embed
      !
      CALL allocate_extforces(embed)
      !
   END SUBROUTINE
   !
   SUBROUTINE allocate_extforces(embed)
      USE kinds,                ONLY : DP
      USE fft_base,             ONLY : dfftp
      USE ions_base,            ONLY : nat
      !
      IMPLICIT NONE
      TYPE(embed_base), INTENT(INOUT) :: embed
      !
      IF (ALLOCATED(embed%extforces)) THEN
         IF (SIZE(embed%extforces,2) /= nat) DEALLOCATE(embed%extforces)
      ENDIF
      IF (.NOT.ALLOCATED(embed%extforces)) THEN
         ALLOCATE(embed%extforces(3, nat))
         embed%extforces = 0.0_DP
      ENDIF
   END SUBROUTINE
   !
   SUBROUTINE arr2pointer_real_1(arr, p, n1)
      USE kinds,                ONLY : DP
      !
      IMPLICIT NONE
      REAL(DP), INTENT(IN), TARGET :: arr(n1)
      REAL(DP), POINTER :: p(:)
      INTEGER,INTENT(IN) :: n1
      !
      p => arr
      !
   END SUBROUTINE
   !
   SUBROUTINE arr2pointer_real_2(arr, p, n1, n2)
      USE kinds,                ONLY : DP
      !
      IMPLICIT NONE
      REAL(DP), INTENT(IN), TARGET :: arr(n1, n2)
      REAL(DP), POINTER :: p(:,:)
      INTEGER,INTENT(IN) :: n1,n2
      !
      p => arr
      !
   END SUBROUTINE
   !
   SUBROUTINE arr2pointer_real_3(arr, p, n1, n2, n3)
      USE kinds,                ONLY : DP
      !
      IMPLICIT NONE
      REAL(DP), INTENT(IN), TARGET :: arr(n1, n2, n3)
      REAL(DP), POINTER :: p(:,:,:)
      INTEGER,INTENT(IN) :: n1,n2,n3
      !
      p => arr
      !
   END SUBROUTINE
   !
   SUBROUTINE arr2pointer_real_4(arr, p, n1, n2, n3, n4)
      USE kinds,                ONLY : DP
      !
      IMPLICIT NONE
      REAL(DP), INTENT(IN), TARGET :: arr(n1, n2, n3, n4)
      REAL(DP), POINTER :: p(:,:,:,:)
      INTEGER,INTENT(IN) :: n1,n2,n3,n4
      !
      p => arr
      !
   END SUBROUTINE
   !
END MODULE qepy_common
