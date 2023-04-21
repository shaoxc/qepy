MODULE qepy_common
   USE kinds,               ONLY : DP
   !
   USE fft_types,           ONLY : fft_type_descriptor
   USE scf,                 ONLY : scf_type
   !
   IMPLICIT NONE
   PRIVATE
   !
   PUBLIC :: arr2pointer, set_embed
   !
#if defined(__MPI)
   logical, public :: is_mpi = .TRUE.
#else
   logical, public :: is_mpi = .FALSE.
#endif
#if defined(_OPENMP)
   logical, public :: is_openmp = .TRUE.
#else
   logical, public :: is_openmp = .FALSE.
#endif
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
   CONTAINS
      !--------------------------------------------------------------------------------
      PROCEDURE :: free => free_tddft
   end type tddft_base
   !
   type, public :: energies_base
      real(kind=dp)                   :: etot           !'Total energy'
      real(kind=dp)                   :: ek             !'Kinetic energy'
      real(kind=dp)                   :: eloc           !'Local energy'
      real(kind=dp)                   :: enl            !'Non-Local energy'
      real(kind=dp)                   :: ewld           !'Ewald energy'
      real(kind=dp)                   :: exc            !'xc contribution'
      real(kind=dp)                   :: ehart          !'hartree energy'
      real(kind=dp)                   :: fock2          !'EXX energy'
      real(kind=dp)                   :: demet          !'Smearing (-TS)'
      real(kind=dp)                   :: elondon        !'Dispersion Correction'
      real(kind=dp)                   :: edftd3         !'DFT-D3 Dispersion'
      real(kind=dp)                   :: exdm           !'Dispersion XDM Correction'
      real(kind=dp)                   :: etsvdw         !'Dispersion T-S Correction'
      real(kind=dp)                   :: eext           !'External forces energy'
      real(kind=dp)                   :: etotefield     !'electric field correction'
      real(kind=dp)                   :: etotgatefield  !'gate field correction'
      real(kind=dp)                   :: eth            !'Hubbard energy'
      real(kind=dp)                   :: epaw           !'one-center paw contrib.'
      real(kind=dp)                   :: ept            !'potentiostat contribution'
      !
      real(kind=dp)                   :: extene         !'External energy0'
      real(kind=dp)                   :: ehf            !'Harris-Foulkes estimate'
      ! some energies details -->
      real(kind=dp)                   :: etxc           ! the exchange and correlation energy
      real(kind=dp)                   :: etxcc          ! the nlcc exchange and correlation
      real(kind=dp)                   :: paw_ehart_ae   !'PAW hartree energy AE'
      real(kind=dp)                   :: paw_ehart_ps   !'PAW hartree energy PS'
      real(kind=dp)                   :: paw_exc_ae     !'PAW xc energy AE'
      real(kind=dp)                   :: paw_exc_ps     !'PAW xc energy PS'
      ! <--
   CONTAINS
      !--------------------------------------------------------------------------------
      PROCEDURE :: reset => assignment_energies
      GENERIC :: assignment(=) => reset
   end type energies_base
   !
   type, public :: embed_base
      type(input_base)                :: input
      type(tddft_base)                :: tddft
      type(energies_base)             :: energies
      real(kind=dp), allocatable      :: extpot(:,:)
      real(kind=dp)                   :: extene = 0.d0
      integer                         :: exttype = 0
      real(kind=dp), allocatable      :: extforces(:,:)
      real(kind=dp)                   :: extstress(3,3) = 0.d0
      logical                         :: initial = .true.
      real(kind=dp)                   :: mix_coef = -1.d0
      logical                         :: finish = .false.
      real(kind=dp)                   :: etotal = 0.d0
      real(kind=dp)                   :: dnorm = 1.d0
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
      !!
      !!
      type(fft_type_descriptor),pointer :: dfftp
      type(fft_type_descriptor),pointer :: dffts
      type(scf_type),pointer            :: rho
      type(scf_type),pointer            :: v
      type(scf_type),pointer            :: vnew
   CONTAINS
      !--------------------------------------------------------------------------------
      PROCEDURE :: allocate_extpot
      PROCEDURE :: allocate_extforces
      PROCEDURE :: free => free_embed
   end type embed_base
   !
   !
   TYPE ( embed_base ), public, pointer :: embed
   TYPE ( embed_base ), public, target  :: messenger
   !
   !
   INTERFACE arr2pointer
      MODULE PROCEDURE arr2pointer_real_1, arr2pointer_real_2, arr2pointer_real_3, arr2pointer_real_4
   END INTERFACE
   !
CONTAINS
   !
   SUBROUTINE set_embed(obj)
      !
      USE fft_base,             ONLY : dfftp, dffts
      USE scf,                  ONLY : rho, vnew, v
      IMPLICIT NONE
      TYPE(embed_base), INTENT(INOUT), TARGET :: obj
      !
      embed => obj
      call set_target_dfft(obj, dfftp, dffts)
      call set_target_scf(obj, rho, vnew, v)
      !
   END SUBROUTINE
   !
   SUBROUTINE set_target_dfft(obj, dfftp, dffts)
      USE fft_types,           ONLY : fft_type_descriptor
      !
      IMPLICIT NONE
      TYPE(fft_type_descriptor), INTENT(IN), TARGET :: dfftp, dffts
      TYPE(embed_base), INTENT(INOUT) :: obj
      !
      obj%dfftp => dfftp
      obj%dffts => dffts
      !
   END SUBROUTINE
   !
   SUBROUTINE set_target_scf(obj, rho, v, vnew)
      USE scf,                 ONLY : scf_type
      !
      IMPLICIT NONE
      TYPE(scf_type), INTENT(IN), TARGET :: rho, v, vnew
      TYPE(embed_base), INTENT(INOUT) :: obj
      !
      obj%rho => rho
      obj%v => v
      obj%vnew => vnew
      !
   END SUBROUTINE
   !
   SUBROUTINE free_embed(obj)
      !
      IMPLICIT NONE
      CLASS(embed_base), INTENT(INOUT) :: obj
      !
      IF (ALLOCATED(obj%extpot)) DEALLOCATE(obj%extpot)
      IF (ALLOCATED(obj%extforces)) DEALLOCATE(obj%extforces)
      call obj%tddft%free()
      !
   END SUBROUTINE
   !
   SUBROUTINE free_tddft(obj)
      !
      IMPLICIT NONE
      CLASS(tddft_base), INTENT(INOUT) :: obj
      !
      IF (ALLOCATED(obj%dipole)) DEALLOCATE(obj%dipole)
      !
   END SUBROUTINE
   !
   SUBROUTINE allocate_extpot(embed)
      USE kinds,                ONLY : DP
      USE fft_base,             ONLY : dfftp
      USE lsda_mod,             ONLY : lsda, nspin
      !
      IMPLICIT NONE
      CLASS(embed_base), INTENT(INOUT) :: embed
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
   SUBROUTINE allocate_extforces(embed)
      USE kinds,                ONLY : DP
      USE fft_base,             ONLY : dfftp
      USE ions_base,            ONLY : nat
      !
      IMPLICIT NONE
      CLASS(embed_base), INTENT(INOUT) :: embed
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
   SUBROUTINE assignment_energies(obj, value)
      !
      IMPLICIT NONE
      CLASS(energies_base),INTENT(INOUT) :: obj
      REAL(DP),INTENT(IN)                :: value
       
      obj%etot          = value
      obj%ek            = value
      obj%eloc          = value
      obj%enl           = value
      obj%ewld          = value
      obj%exc           = value
      obj%ehart         = value
      obj%fock2         = value
      obj%demet         = value
      obj%elondon       = value
      obj%edftd3        = value
      obj%exdm          = value
      obj%etsvdw        = value
      obj%eext          = value
      obj%etotefield    = value
      obj%etotgatefield = value
      obj%eth           = value
      obj%epaw          = value
      obj%ept           = value
      obj%extene        = value
      obj%ehf           = value
      obj%etxc          = value
      obj%etxcc         = value
      obj%paw_ehart_ae  = value
      obj%paw_ehart_ps  = value
      obj%paw_exc_ae    = value
      obj%paw_exc_ps    = value
      
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
