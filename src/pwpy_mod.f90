MODULE pwpy_mod
   USE kinds,                ONLY : DP
   IMPLICIT NONE
   PUBLIC
   REAL(DP), ALLOCATABLE :: scf__vnew__of_r(:,:)
   CONTAINS
      SUBROUTINE pwpy_init_pointer()
         use scf, only: rho,v,vnew
         IMPLICIT NONE
         REAL(DP), POINTER :: scf__rho__of_r(:,:)
         REAL(DP), POINTER :: scf__v__of_r(:,:)
         REAL(DP), POINTER :: scf__vnew__of_r(:,:)
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
         IMPLICIT NONE
         !REAL(DP), INTENT(OUT) :: rhor(dfftp%nnr,nspin)
         REAL(DP), INTENT(OUT) :: rhor(:,:)
         rhor(:,:)=rho%of_r(:,:)
      END SUBROUTINE 

      SUBROUTINE pwpy_set_rho(rhor)
         USE kinds,                ONLY : DP
         USE fft_rho,              ONLY : rho_g2r, rho_r2g
         USE fft_base,         ONLY : dfftp, dffts
         use scf, only: rho !! the charge density and its other components
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: rhor(:,:)
         rho%of_r(:,:)=rhor
         CALL rho_r2g(dfftp, rho%of_r, rho%of_g )
      END SUBROUTINE 

      SUBROUTINE pwpy_set_v(vin)
         USE kinds,                ONLY : DP
         USE fft_rho,              ONLY : rho_g2r, rho_r2g
         USE fft_base,         ONLY : dfftp, dffts
         use scf, only: rho, &!! the charge density and its other components
            v,         & !! the scf potential
            vnew,      & !! used to correct the forces
            v_of_0,    & !! vltot(G=0)
            vltot,     & !! the local potential in real space
            vrs,       & !! the total pot. in real space (smooth grid)
            rho_core,  & !! the core charge in real space
            kedtau,    & !! position dependent kinetic energy enhancement factor
            rhog_core !! the core charge in reciprocal space

         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: vin(:,:)
      END SUBROUTINE 

      SUBROUTINE pwpy_get_grid(nr)
         USE kinds,                ONLY : DP
         USE fft_base,         ONLY : dfftp, dffts
         IMPLICIT NONE
         INTEGER, INTENT(OUT) :: nr(3)
         nr=(/dfftp%nr1, dfftp%nr2, dfftp%nr3/)
      END SUBROUTINE 

      SUBROUTINE pwpy_get_forces(forces)
         USE kinds,                ONLY : DP
         USE force_mod,         ONLY : force, lforce, sumfor
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: forces(:,:)
      END SUBROUTINE 
END MODULE pwpy_mod
