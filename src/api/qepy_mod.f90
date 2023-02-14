MODULE qepy_mod
   USE kinds,                   ONLY : DP
   USE qepy_scatter_mod,        ONLY : gather_grid, scatter_grid
   USE qepy_common,             ONLY : embed_base, input_base
   USE fft_base,                ONLY : dfftp
   !
   IMPLICIT NONE
   PUBLIC
   !
   INTERFACE mp_gather
      MODULE PROCEDURE mp_gather_real, mp_gather_complex
   END INTERFACE
   !
   INTERFACE mp_scatter
      MODULE PROCEDURE mp_scatter_real, mp_scatter_complex
   END INTERFACE
   !
   INTERFACE qepy_get_value
      MODULE PROCEDURE qepy_get_value_real_1, qepy_get_value_real_2
   END INTERFACE
   !
CONTAINS

   SUBROUTINE mp_gather_real(fin, fout)
      USE kinds,                ONLY : DP
      USE fft_base,             ONLY : dfftp
      !
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: fin(:)
      REAL(DP), INTENT(OUT) :: fout(:)
      !
      IF (dfftp%nproc > 1) THEN
         CALL gather_grid(dfftp, fin, fout)
      ELSE
         fout(:) = fin(:)
      ENDIF
   END SUBROUTINE

   SUBROUTINE mp_scatter_real(fin, fout)
      USE kinds,                ONLY : DP
      USE fft_base,             ONLY : dfftp
      !
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: fin(:)
      REAL(DP), INTENT(OUT) :: fout(:)
      !
      IF (dfftp%nproc > 1) THEN
         CALL scatter_grid(dfftp, fin, fout)
      ELSE
         fout(:) = fin(:)
      ENDIF
   END SUBROUTINE

   SUBROUTINE mp_gather_complex(fin, fout)
      USE kinds,                ONLY : DP
      USE fft_base,             ONLY : dfftp
      !
      IMPLICIT NONE
      COMPLEX(DP), INTENT(IN) :: fin(:)
      COMPLEX(DP), INTENT(OUT) :: fout(:)
      !
      IF (dfftp%nproc > 1) THEN
         CALL gather_grid(dfftp, fin, fout)
      ELSE
         fout(:) = fin(:)
      ENDIF
   END SUBROUTINE

   SUBROUTINE mp_scatter_complex(fin, fout)
      USE kinds,                ONLY : DP
      USE fft_base,             ONLY : dfftp
      !
      IMPLICIT NONE
      COMPLEX(DP), INTENT(IN) :: fin(:)
      COMPLEX(DP), INTENT(OUT) :: fout(:)
      !
      IF (dfftp%nproc > 1) THEN
         CALL scatter_grid(dfftp, fin, fout)
      ELSE
         fout(:) = fin(:)
      ENDIF
   END SUBROUTINE

   SUBROUTINE qepy_get_value_real_1(fin, fout, gather, scatter)
      USE kinds,                ONLY : DP
      USE fft_base,             ONLY : dfftp
      !
      IMPLICIT NONE
      real(DP), INTENT(IN)        :: fin(:)
      real(DP), INTENT(OUT)       :: fout(:)
      LOGICAL,INTENT(in),OPTIONAL :: gather
      LOGICAL,INTENT(in),OPTIONAL :: scatter
      !
      INTEGER                     :: nnr
      LOGICAL                     :: gather_ = .FALSE.
      LOGICAL                     :: scatter_ = .FALSE.
      !
      IF ( present(gather) ) gather_ = gather
      IF ( present(scatter) ) scatter_ = scatter
      !
      IF ( gather ) THEN
         CALL mp_gather(fin, fout)
      ELSE IF ( scatter ) THEN
         CALL mp_scatter(fin, fout)
      ELSE
         !nnr = dfftp%nr1x* dfftp%my_nr2p* dfftp%my_nr3p
         nnr = size(fin)
         fout(1:nnr) = fin
         fout(nnr:size(fout)) = 0.0_DP
      ENDIF
      !
   END SUBROUTINE

   SUBROUTINE qepy_get_value_real_2(fin, fout, gather, scatter)
      USE kinds,                ONLY : DP
      !
      IMPLICIT NONE
      real(DP), INTENT(IN)        :: fin(:,:)
      real(DP), INTENT(OUT)       :: fout(:,:)
      LOGICAL,INTENT(in),OPTIONAL :: gather
      LOGICAL,INTENT(in),OPTIONAL :: scatter
      !
      INTEGER                     :: ispin, nspin
      LOGICAL                     :: gather_ = .FALSE.
      LOGICAL                     :: scatter_ = .FALSE.
      !
      IF ( present(gather) ) gather_ = gather
      IF ( present(scatter) ) scatter_ = scatter
      !
      nspin = size(fin, 2)
      DO ispin = 1, nspin
         call qepy_get_value_real_1(fin(:, ispin), fout(:, ispin), gather_, scatter_)
      END DO
      !
   END SUBROUTINE

   SUBROUTINE qepy_get_rho(rhor, gather)
      USE kinds,                ONLY : DP
      use scf,                  ONLY : rho, rhoz_or_updw
      USE fft_base,             ONLY : dfftp, dffts
      USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
      !
      IMPLICIT NONE
      REAL(DP), INTENT(OUT) :: rhor(:,:)
      LOGICAL,INTENT(in),OPTIONAL :: gather
      !
      LOGICAL :: gather_ = .true.
      !
      IF ( present(gather) ) gather_ = gather
      !
      IF (nspin > 1) CALL rhoz_or_updw( rho, 'only_r', '->updw' )
      !
      call qepy_get_value(rho%of_r, rhor, gather = gather_)
      !
      IF (nspin > 1) CALL rhoz_or_updw( rho, 'only_r', '->rhoz' )
   END SUBROUTINE

   SUBROUTINE qepy_set_rho(rhor, gather)
      USE kinds,                ONLY : DP
      USE fft_rho,              ONLY : rho_g2r, rho_r2g
      USE fft_base,             ONLY : dfftp, dffts
      use scf,                  ONLY : rho, rhoz_or_updw
      USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
      !
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: rhor(:,:)
      LOGICAL,INTENT(in),OPTIONAL :: gather
      !
      LOGICAL :: gather_ = .true.
      !
      IF ( present(gather) ) gather_ = gather
      !
      IF (nspin > 1) CALL rhoz_or_updw( rho, 'only_r', '->updw' )
      !
      call qepy_get_value(rhor, rho%of_r, scatter = gather_)
      !
      IF (nspin > 1) CALL rhoz_or_updw( rho, 'only_r', '->rhoz' )
      CALL rho_r2g(dfftp, rho%of_r, rho%of_g )
   END SUBROUTINE

   SUBROUTINE qepy_get_rho_core(rhoc, gather)
      USE kinds,                ONLY : DP
      use scf,                  ONLY : rho_core !! the core charge in real space
      USE fft_base,             ONLY : dfftp, dffts
      IMPLICIT NONE
      REAL(DP), INTENT(OUT) :: rhoc(:)
      LOGICAL,INTENT(in),OPTIONAL :: gather
      !
      LOGICAL :: gather_ = .true.
      !
      IF ( present(gather) ) gather_ = gather
      !
      call qepy_get_value(rho_core, rhoc, gather = gather_)
      !
   END SUBROUTINE

   SUBROUTINE qepy_set_rho_core(rhoc, gather)
      USE kinds,                ONLY : DP
      use scf,                  ONLY : rho_core !! the core charge in real space
      USE fft_base,             ONLY : dfftp, dffts
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: rhoc(:)
      LOGICAL,INTENT(in),OPTIONAL :: gather
      !
      LOGICAL :: gather_ = .true.
      !
      IF ( present(gather) ) gather_ = gather
      !
      call qepy_get_value(rhoc, rho_core, scatter = gather_)
      !
   END SUBROUTINE

   SUBROUTINE qepy_set_extpot(embed, vin, gather)
      USE kinds,                ONLY : DP
      USE fft_rho,              ONLY : rho_g2r, rho_r2g
      USE fft_base,             ONLY : dfftp, dffts
      USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
      USE mp,                   ONLY : mp_bcast
      !
      IMPLICIT NONE
      TYPE(embed_base), INTENT(INOUT) :: embed
      REAL(DP), INTENT(IN) :: vin(:,:)
      LOGICAL,INTENT(in),OPTIONAL :: gather
      !
      INTEGER :: ispin, ns
      LOGICAL :: gather_ = .true.
      !
      IF ( present(gather) ) gather_ = gather
      !
      call embed%allocate_extpot()
      !
      ns = size(vin,2)
      CALL mp_bcast(ns, dfftp%root, dfftp%comm)
      !
      call qepy_get_value(vin(:,1:ns), embed%extpot(:,1:ns), scatter = gather_)
      !
      DO ispin = ns+1, nspin
         embed%extpot(:,ispin) = embed%extpot(:,1)
      END DO
      !
   END SUBROUTINE

   FUNCTION qepy_get_grid(nr, gather) RESULT( nrw )
      USE kinds,                ONLY : DP
      USE fft_base,             ONLY : dfftp
      !
      IMPLICIT NONE
      INTEGER,INTENT(OUT),OPTIONAL :: nr(3)
      LOGICAL,INTENT(in),OPTIONAL  :: gather
      !
      INTEGER                      :: nrw(3)
      !
      IF ( present(gather) ) THEN
         nrw = qepy_get_grid_shape(dfftp, gather)
      ELSE
         nrw = qepy_get_grid_shape(dfftp)
      ENDIF
      !
      IF ( present(nr) ) nr = nrw
   END FUNCTION

   FUNCTION qepy_get_grid_shape(dfft, gather) RESULT( nrw )
      USE kinds,                ONLY : DP
      USE fft_types,            ONLY : fft_type_descriptor
      !
      IMPLICIT NONE
      TYPE(fft_type_descriptor),INTENT(IN) :: dfft
      LOGICAL,INTENT(in),OPTIONAL          :: gather
      !
      LOGICAL                              :: gather_ = .true.
      INTEGER                              :: nrw(3)
      !
      IF ( present(gather) ) gather_ = gather
      !
      IF ( gather_ ) THEN
         nrw =(/dfft%nr1, dfft%nr2, dfft%nr3/)
      ELSE
         nrw =(/dfft%nr1x, dfft%my_nr2p, dfft%my_nr3p/)
      ENDIF
      !
   END FUNCTION

   FUNCTION qepy_get_grid_smooth(nr, gather) RESULT(nrw)
      USE kinds,                ONLY : DP
      USE fft_base,             ONLY : dffts
      !
      IMPLICIT NONE
      INTEGER,INTENT(OUT),OPTIONAL :: nr(3)
      LOGICAL,INTENT(in),OPTIONAL  :: gather
      !
      INTEGER                      :: nrw(3)
      !
      IF ( present(gather) ) THEN
         nrw = qepy_get_grid_shape(dffts, gather)
      ELSE
         nrw = qepy_get_grid_shape(dffts)
      ENDIF
      !
      IF ( present(nr) ) nr = nrw
   END FUNCTION

   SUBROUTINE qepy_set_stdout(fname, uni, append)
      USE io_global,            ONLY : stdout, ionode
      !
      INTEGER                  :: ierr
      CHARACTER(LEN=*),INTENT(IN),OPTIONAL  :: fname
      INTEGER,INTENT(in),OPTIONAL :: uni
      LOGICAL,INTENT(in),OPTIONAL :: append
      !
      LOGICAL :: exst
      !
      IF ( .NOT. present(fname) ) return
      IF ( present(uni) ) THEN
         stdout = uni
      ELSE
         stdout = 666
      ENDIF
      !
      exst = .false.
      IF ( present(append)) THEN
         IF (append) INQUIRE (file = TRIM(fname), exist = exst)
      ENDIF
      !
      IF (exst) THEN
         OPEN (UNIT = stdout, FILE = TRIM(fname), FORM = 'formatted', POSITION = 'append', iostat = ierr )
      ELSE
         OPEN (UNIT = stdout, FILE = TRIM(fname), FORM = 'formatted', STATUS = 'unknown', iostat = ierr )
      ENDIF
   END SUBROUTINE

   SUBROUTINE qepy_write_stdout(fstr)
      USE io_global,            ONLY : stdout, ionode
      !
      INTEGER                  :: ierr
      CHARACTER(LEN=*),INTENT(IN)  :: fstr
      !
      IF(ionode) WRITE(stdout,'(A)') fstr
   END SUBROUTINE

   SUBROUTINE qepy_close_stdout(fname)
      USE io_global,            ONLY : stdout, ionode
      !
      INTEGER                  :: ierr
      CHARACTER(LEN=*),INTENT(IN)  :: fname
      !
      IF(ionode) close(stdout)
   END SUBROUTINE

   SUBROUTINE qepy_get_evc(ik, wfc)
      USE kinds,                ONLY : DP
      USE io_files,             ONLY : iunwfc, nwordwfc
      USE buffers,              ONLY : get_buffer
      USE wavefunctions,        ONLY : evc
      USE klist,                ONLY : nks
      !
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: ik
      COMPLEX(DP), INTENT(OUT),OPTIONAL :: wfc(:,:)
      !
      IF ( nks > 1 ) CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
      IF ( present(wfc) ) THEN
         wfc = evc
      ENDIF
   END SUBROUTINE

   SUBROUTINE qepy_get_wf(ik, ibnd, wf, gather)
      USE kinds,                ONLY : DP
      USE io_files,             ONLY : iunwfc, nwordwfc
      USE buffers,              ONLY : get_buffer
      USE wavefunctions,        ONLY : evc, psic
      USE klist,                ONLY : nks, igk_k, ngk
      USE fft_base,             ONLY : dfftp, dffts
      USE fft_interfaces,       ONLY : invfft
      USE control_flags,        ONLY : diago_full_acc, gamma_only, lxdm, tqr
      USE bp,                   ONLY : lelfield
      !
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: ik, ibnd
      COMPLEX(DP), INTENT(OUT) :: wf(:)
      LOGICAL,INTENT(in),OPTIONAL :: gather
      !
      INTEGER :: j, nnr, npw
      LOGICAL :: gather_ = .true.
      !
      IF ( present(gather) ) gather_ = gather
      !
      IF ( dffts%has_task_groups ) THEN
         call errore('qepy_get_wf', 'Sorry this one not support task-group version', 1)
      ENDIF
      !
      IF ( nks > 1 .OR. lelfield ) CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
      !$omp parallel
      psic(:) = (0.0_DP, 0.0_DP)
      npw = ngk(ik)
      !$omp do
      IF ( gamma_only ) THEN
         psic(dffts%nl (1:npw))  = evc(1:npw,ibnd)
         psic(dffts%nlm(1:npw)) = CONJG( evc(1:npw,ibnd) )
      ELSE
         DO j = 1, npw
            psic(dffts%nl(igk_k(j,ik))) = evc(j,ibnd)
         ENDDO
      END IF
      !$omp end do nowait
      !$omp end parallel
      CALL invfft ('Wave', psic, dffts)
      !
      IF ( gather_ ) THEN
         CALL mp_gather(psic(1:dffts%nnr), wf)
      ELSE
         nnr = min(size(wf), dffts%nnr)
         wf(1:nnr) = psic(1:nnr)
         wf(nnr:size(wf)) = (0.0_DP, 0.0_DP)
      ENDIF
   END SUBROUTINE

   SUBROUTINE qepy_get_vkb(ik, vk, gather)
      USE kinds,                ONLY : DP
      USE io_files,             ONLY : iunwfc, nwordwfc
      USE buffers,              ONLY : get_buffer
      USE wavefunctions,        ONLY : evc, psic
      USE klist,                ONLY : nks, igk_k, ngk, xk
      USE uspp,                 ONLY : vkb, nkb
      USE fft_base,             ONLY : dfftp, dffts
      USE fft_interfaces,       ONLY : invfft
      USE control_flags,        ONLY : diago_full_acc, gamma_only, lxdm, tqr
      !
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: ik
      COMPLEX(DP), INTENT(OUT) :: vk(:,:)
      LOGICAL,INTENT(in),OPTIONAL :: gather
      !
      INTEGER :: i, j, nnr, npw
      LOGICAL :: gather_ = .true.
      !
      IF ( present(gather) ) gather_ = gather
      !
      IF ( dffts%has_task_groups ) THEN
         call errore('qepy_get_vkb', 'Sorry this one not support task-group version', 1)
      ENDIF
      !
      IF ( nkb > 0 ) CALL init_us_2( ngk(ik), igk_k(1,ik), xk(1,ik), vkb )
      !
      vk(:,:) = (0.0_DP, 0.0_DP)
      DO i = 1, nkb
         !$omp parallel
         psic(:) = (0.0_DP, 0.0_DP)
         npw = ngk(ik)
         !$omp do
         IF ( gamma_only ) THEN
            psic(dffts%nl (1:npw))  = vkb(1:npw,i)
            psic(dffts%nlm(1:npw)) = CONJG( vkb(1:npw,i) )
         ELSE
            DO j = 1, npw
               psic(dffts%nl(igk_k(j,ik))) = vkb(j,i)
            ENDDO
         END IF
         !$omp end do nowait
         !$omp end parallel
         CALL invfft ('Wave', psic, dffts)
         !
         IF ( gather_ ) THEN
            CALL mp_gather(psic(1:dffts%nnr), vk(:, i))
         ELSE
            nnr = min(size(vk, 1), dffts%nnr)
            vk(1:nnr, i) = psic(1:nnr)
         ENDIF
      ENDDO
   END SUBROUTINE

   SUBROUTINE qepy_set_extforces(embed, forces)
      USE kinds,                ONLY : DP
      USE ions_base,            ONLY : nat, ntyp => nsp
      !
      IMPLICIT NONE
      TYPE(embed_base), INTENT(INOUT) :: embed
      REAL(DP), INTENT(IN) :: forces(:,:)
      !
      call embed%allocate_extforces()
      embed%extforces(:,:) = forces(:,1:nat)
      !
   END SUBROUTINE

   SUBROUTINE qepy_calc_effective_potential(embed, potential, gather)
      USE kinds,                ONLY : DP
      USE ions_base,            ONLY : nat, ntyp => nsp
      USE scf,                  ONLY : rho, rho_core, rhog_core, v, vltot, vrs
      USE ener,                 ONLY : etot, hwf_energy, eband, deband, ehart, &
                                       vtxc, etxc, etxcc, ewld, demet, epaw, &
                                       elondon, edftd3, ef_up, ef_dw
      USE ldaU,                 ONLY : eth
      USE extfield,             ONLY : tefield, etotefield
      USE lsda_mod,             ONLY : nspin
      !
      IMPLICIT NONE
      TYPE(embed_base), INTENT(INOUT) :: embed
      REAL(DP), INTENT(OUT),OPTIONAL  :: potential(:,:)
      LOGICAL,INTENT(in),OPTIONAL     :: gather
      !
      REAL(DP) :: charge
      !
      LOGICAL :: gather_ = .true.
      !
      IF ( present(gather) ) gather_ = gather
      !
      CALL qepy_v_of_rho_all( rho, rho_core, rhog_core, &
         ehart, etxc, vtxc, eth, etotefield, charge, v, embed)
      !
      IF ( present(potential) ) THEN
         call qepy_get_value(vrs, potential, gather = gather_)
      ENDIF
      !
   END SUBROUTINE

   SUBROUTINE qepy_set_effective_potential(embed, potential, gather)
      USE kinds,                ONLY : DP
      USE scf,                  ONLY : kedtau, v, vltot, vrs
      USE lsda_mod,             ONLY : nspin
      USE gvecs,                ONLY : doublegrid
      !
      IMPLICIT NONE
      TYPE(embed_base), INTENT(INOUT) :: embed
      REAL(DP),INTENT(IN)             :: potential(:,:)
      LOGICAL,INTENT(in),OPTIONAL     :: gather
      !
      LOGICAL :: gather_ = .true.
      !
      IF ( present(gather) ) gather_ = gather
      !
      call qepy_get_value(potential, vrs, gather = gather_)
      !
      CALL interpolate_vrs( dfftp%nnr, nspin, doublegrid, kedtau, v%kin_r, vrs )
      !
   END SUBROUTINE

   SUBROUTINE qepy_calc_density(rhor, gather)
      USE kinds,                ONLY : DP
      USE wvfct,                ONLY : nbnd, et
      USE klist,                ONLY : nks, nkstot
      !
      IMPLICIT NONE
      REAL(DP), INTENT(OUT),OPTIONAL :: rhor(:,:)
      LOGICAL,INTENT(in),OPTIONAL :: gather
      !
      CALL poolrecover( et, nbnd, nkstot, nks )
      CALL sum_band()
      !
      IF ( present(rhor) ) THEN
         IF ( present(gather) ) THEN
            call qepy_get_rho(rhor, gather)
         ELSE
            call qepy_get_rho(rhor)
         ENDIF
      ENDIF
   END SUBROUTINE

   SUBROUTINE qepy_diagonalize(iter, threshold)
      USE klist,                ONLY : nks, nkstot
      USE bp,                   ONLY : lelfield
      USE control_flags,        ONLY : ethr
      !
      IMPLICIT NONE
      INTEGER,INTENT(in),OPTIONAL :: iter
      REAL(DP),INTENT(in),OPTIONAL :: threshold
      !
      INTEGER                     :: it = 1
      !
      IF ( present(iter) ) it = iter
      IF ( present(threshold) ) ethr = threshold
      !
      IF ( lelfield ) THEN
         CALL c_bands_efield( it )
      ELSE
         CALL c_bands( it )
      ENDIF
   END SUBROUTINE

END MODULE qepy_mod
