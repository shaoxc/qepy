MODULE pwpy_embed
   !
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
   end type embed_base
   !
END MODULE pwpy_embed
