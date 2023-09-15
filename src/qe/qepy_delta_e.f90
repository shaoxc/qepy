!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!qepy copy from electrons.f90
     FUNCTION qepy_delta_e(vr)
       !-----------------------------------------------------------------------
       ! This function computes delta_e, where:
       !
       ! ... delta_e =  - \int rho%of_r(r)  v%of_r(r)
       !                - \int rho%kin_r(r) v%kin_r(r) [for Meta-GGA]
       !                - \sum rho%ns       v%ns       [for DFT+Hubbard]
       !                - \sum becsum       D1_Hxc     [for PAW]
       !
       USE xc_lib,  ONLY : xclib_dft_is
       !
  USE kinds,                ONLY : DP
  USE check_stop,           ONLY : check_stop_now, stopped_by_user
  USE io_global,            ONLY : stdout, ionode
  USE cell_base,            ONLY : at, bg, alat, omega, tpiba2
  USE ions_base,            ONLY : zv, nat, nsp, ityp, tau, compute_eextfor, atm, &
                                   ntyp => nsp
  USE basis,                ONLY : starting_pot
  USE bp,                   ONLY : lelfield
  USE fft_base,             ONLY : dfftp
  USE gvect,                ONLY : ngm, gstart, g, gg, gcutm
  USE gvecs,                ONLY : doublegrid, ngms
  USE klist,                ONLY : xk, wk, nelec, ngk, nks, nkstot, lgauss, &
                                   two_fermi_energies, tot_charge
  USE fixed_occ,            ONLY : one_atom_occupations
  USE lsda_mod,             ONLY : lsda, nspin, magtot, absmag, isk
  USE vlocal,               ONLY : strf
  USE wvfct,                ONLY : nbnd, et
  USE gvecw,                ONLY : ecutwfc
  USE ener,                 ONLY : etot, hwf_energy, eband, deband, ehart, &
                                   vtxc, etxc, etxcc, ewld, demet, epaw, &
                                   elondon, edftd3, ef_up, ef_dw, exdm, ef, &
                                   egrand, vsol, esol, esic, esci
  USE scf,                  ONLY : scf_type, scf_type_COPY, bcast_scf_type,&
                                   create_scf_type, destroy_scf_type, &
                                   open_mix_file, close_mix_file, &
                                   rho, rho_core, rhog_core, v, vltot, vrs, &
                                   kedtau, vnew
  USE control_flags,        ONLY : mixing_beta, tr2, ethr, niter, nmix, &
                                   iprint, conv_elec, sic, &
                                   restart, io_level, do_makov_payne,  &
                                   gamma_only, iverbosity, textfor,     &
                                   llondon, ldftd3, scf_must_converge, lxdm, ts_vdw, &
                                   mbd_vdw, use_gpu
  USE control_flags,        ONLY : n_scf_steps, scf_error, scissor
  USE sci_mod,              ONLY : sci_iter

  USE io_files,             ONLY : iunmix, output_drho
  USE ldaU,                 ONLY : eth, lda_plus_u, lda_plus_u_kind, &
                                   niter_with_fixed_ns, hub_pot_fix, &
                                   nsg, nsgnew, v_nsg, at_sc, neighood, &
                                   ldim_u, is_hubbard_back
  USE extfield,             ONLY : tefield, etotefield, gate, etotgatefield !TB
  USE noncollin_module,     ONLY : noncolin, magtot_nc, i_cons,  bfield, &
                                   lambda, report, domag, nspin_mag
  USE io_rho_xml,           ONLY : write_scf
  USE uspp,                 ONLY : okvan
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp_pools,             ONLY : root_pool, me_pool, my_pool_id, &
                                   inter_pool_comm, intra_pool_comm
  USE mp,                   ONLY : mp_sum, mp_bcast
  !
  USE london_module,        ONLY : energy_london
  USE dftd3_api,            ONLY : dftd3_pbc_dispersion, get_atomic_number
  USE dftd3_qe,             ONLY : dftd3
  USE xdm_module,           ONLY : energy_xdm
  USE tsvdw_module,         ONLY : EtsvdW
  !
  USE paw_variables,        ONLY : okpaw, ddd_paw, total_core_energy, only_paw
       !
       IMPLICIT NONE
       !
       REAL(DP) :: delta_e
       REAL(DP) :: delta_e_hub
       INTEGER  :: ir, na1, nt1, na2, nt2, m1, m2, equiv_na2, viz, is
       !
       REAL(DP) :: qepy_delta_e
       REAL(DP) :: vr(size(vrs,1),size(vrs,2))
       LOGICAL  :: lhb
       INTEGER  :: nt
       !
       lhb = .FALSE.
       IF ( lda_plus_u )  THEN
          DO nt = 1, ntyp
             IF (is_hubbard_back(nt)) lhb = .TRUE.
          ENDDO
       ENDIF
       !
       delta_e = 0._DP
       IF ( nspin==2 ) THEN
          !
          DO ir = 1,dfftp%nnr
            delta_e = delta_e - ( rho%of_r(ir,1) + rho%of_r(ir,2) ) * vr(ir,1) &  ! up
                              - ( rho%of_r(ir,1) - rho%of_r(ir,2) ) * vr(ir,2)    ! dw
          ENDDO 
          delta_e = 0.5_DP*delta_e
          !
       ELSE
          delta_e = - SUM( rho%of_r(:,1:nspin_mag)*vr(:,1:nspin_mag) )
       ENDIF
       !
       IF ( xclib_dft_is('meta') ) &
          delta_e = delta_e - SUM( rho%kin_r(:,1:nspin_mag)*v%kin_r(:,1:nspin_mag) )
       !
       delta_e = omega * delta_e / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
       !
       CALL mp_sum( delta_e, intra_bgrp_comm )
       !
       IF (lda_plus_u .AND. (.NOT.hub_pot_fix)) THEN
         IF (lda_plus_u_kind.EQ.0) THEN
            delta_e_hub = - SUM( rho%ns(:,:,:,:)*v%ns(:,:,:,:) )
            IF (lhb) delta_e_hub = delta_e_hub - SUM(rho%nsb(:,:,:,:)*v%nsb(:,:,:,:))
            IF (nspin==1) delta_e_hub = 2.d0 * delta_e_hub
            delta_e = delta_e + delta_e_hub
         ELSEIF (lda_plus_u_kind.EQ.1) THEN
            IF (noncolin) THEN
              delta_e_hub = - SUM( rho%ns_nc(:,:,:,:)*v%ns_nc(:,:,:,:) )
              delta_e = delta_e + delta_e_hub
            ELSE
              delta_e_hub = - SUM( rho%ns(:,:,:,:)*v%ns(:,:,:,:) )
              IF (nspin==1) delta_e_hub = 2.d0 * delta_e_hub
              delta_e = delta_e + delta_e_hub
            ENDIF
         ELSEIF (lda_plus_u_kind.EQ.2) THEN
            delta_e_hub = 0.0_dp
            DO is = 1, nspin
               DO na1 = 1, nat
                  nt1 = ityp(na1)
                  DO viz = 1, neighood(na1)%num_neigh
                     na2 = neighood(na1)%neigh(viz)
                     equiv_na2 = at_sc(na2)%at
                     nt2 = ityp(equiv_na2)
                     IF ( ANY(v_nsg(:,:,viz,na1,is).NE.0.0d0) ) THEN
                        DO m1 = 1, ldim_u(nt1)
                           DO m2 = 1, ldim_u(nt2)
                              delta_e_hub = delta_e_hub - &
                                  nsgnew(m2,m1,viz,na1,is)*v_nsg(m2,m1,viz,na1,is)
                           ENDDO
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
            IF (nspin==1) delta_e_hub = 2.d0 * delta_e_hub
            delta_e = delta_e + delta_e_hub
         ENDIF
       ENDIF
       !
       IF (okpaw) delta_e = delta_e - SUM( ddd_paw(:,:,:)*rho%bec(:,:,:) )
       !
       qepy_delta_e = delta_e
       !
       RETURN
       !
     END FUNCTION qepy_delta_e
