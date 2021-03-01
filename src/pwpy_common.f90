MODULE pwpy_common
   USE kinds,                ONLY : DP
   IMPLICIT NONE
   PUBLIC
   !
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
   !
   type, public :: input_base
      INTEGER            :: my_world_comm = 0
      LOGICAL            :: start_images = .false.
      CHARACTER(len=256) :: filename = ''
      CHARACTER(len=256) :: code = 'PWPY'
      CHARACTER(len=256) :: tmp_dir = './'
      CHARACTER(len=256) :: wfc_dir = 'undefined'
      CHARACTER(len=256) :: prefix  = 'os'
      ! ...  If needwf=.t. performs wavefunction-related initialization as well
      LOGICAL            :: needwf
   end type input_base

END MODULE pwpy_common
