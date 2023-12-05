!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE qepy_v_of_rho_all( rho, rho_core, rhog_core, &
                     ehart, etxc, vtxc, eth, etotefield, charge, v )
  !----------------------------------------------------------------------------
  !! This routine computes the Hartree and Exchange and Correlation
  !! potential and energies which corresponds to a given charge density
  !! The XC potential is computed in real space, while the
  !! Hartree potential is computed in reciprocal space.
  !
  USE kinds,            ONLY : DP
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : ngm
  USE noncollin_module, ONLY : noncolin, nspin_lsda
  USE ions_base,        ONLY : nat, tau
  USE ldaU,             ONLY : lda_plus_u, lda_plus_u_kind, ldmx_b, &
                               nsg, v_nsg 
  USE xc_lib,           ONLY : xclib_dft_is
  USE scf,              ONLY : scf_type
  USE cell_base,        ONLY : alat
  USE io_global,        ONLY : stdout
  USE control_flags,    ONLY : ts_vdw, mbd_vdw, sic
  USE tsvdw_module,     ONLY : tsvdw_calculate, UtsvdW
  USE libmbd_interface, ONLY : mbd_interface
  USE sic_mod,          ONLY : add_vsic
  !
  USE scf,                  ONLY : vltot, vrs, kedtau, vnew
  USE gvecs,                ONLY : doublegrid
  USE lsda_mod,             ONLY : lsda, nspin, magtot, absmag, isk
  USE plugin_variables,     ONLY : plugin_etot
  USE dfunct,               ONLY : newd
  USE paw_variables,        ONLY : okpaw, ddd_paw, total_core_energy, only_paw
  USE paw_onecenter,        ONLY : PAW_potential
  USE paw_symmetry,         ONLY : PAW_symmetrize_ddd
  USE ener,                 ONLY : epaw
  USE control_flags,        ONLY : use_gpu
  USE scf_gpum,             ONLY : using_vrs
  USE dfunct_gpum,          ONLY : newd_gpu
  !
  USE qepy_common,      ONLY : embed
  !
  IMPLICIT NONE
  !
  TYPE(scf_type), INTENT(INOUT) :: rho
  !! the valence charge
  TYPE(scf_type), INTENT(INOUT) :: v
  !! the scf (Hxc) potential 
  !=================> NB: NOTE that in F90 derived data type must be INOUT and 
  !=================> not just OUT because otherwise their allocatable or pointer
  !=================> components are NOT defined 
  REAL(DP), INTENT(IN) :: rho_core(dfftp%nnr)
  !! the core charge
  COMPLEX(DP), INTENT(IN) :: rhog_core(ngm)
  !! the core charge in reciprocal space
  REAL(DP), INTENT(OUT) :: vtxc
  !! the integral V_xc * rho
  REAL(DP), INTENT(OUT) :: etxc
  !! the E_xc energy
  REAL(DP), INTENT(OUT) :: ehart
  !! the hartree energy
  REAL(DP), INTENT(OUT) :: eth
  !! the hubbard energy
  REAL(DP), INTENT(OUT) :: charge
  !! the integral of the charge
  REAL(DP) :: eth1
  !! the hubbard energy coming from the background states
  REAL(DP), INTENT(INOUT) :: etotefield
  !! electric field energy - inout due to the screwed logic of add_efield
  !
  REAL(DP) :: etot_cmp_paw(nat,2,2)
  !
  call qepy_v_of_rho( rho, rho_core, rhog_core, &
                     ehart, etxc, vtxc, eth, etotefield, charge, v)
  IF (okpaw) THEN
     CALL PAW_potential( rho%bec, ddd_paw, epaw, etot_cmp_paw )
     CALL PAW_symmetrize_ddd( ddd_paw )
  ENDIF
     !
!#if defined (__LEGACY_PLUGINS) 
     !CALL plugin_scf_energy (plugin_etot, rhoin) 
     !!
     !CALL plugin_scf_potential(rhoin, conv_elec, dr2, vltot)
!#endif 
!#if defined (__ENVIRON)
     !IF (use_environ) THEN
        !CALL calc_environ_energy(plugin_etot, .TRUE.)
        !CALL calc_environ_potential(rhoin, conv_elec, dr2, vltot)
     !END IF
!#endif
!#if defined (__OSCDFT)
     !IF (use_oscdft) THEN
        !CALL oscdft_scf_energy(oscdft_ctx, plugin_etot)
     !END IF
!#endif
     !
     ! ... define the total local potential (external + scf)
     !
     IF (ALLOCATED(embed%extpot)) v%of_r = v%of_r + embed%extpot
     !
     CALL using_vrs(1)
     CALL sum_vrs( dfftp%nnr, nspin, vltot, v%of_r, vrs )
     !
     ! ... interpolate the total local potential
     !
     CALL using_vrs(1) ! redundant
     CALL interpolate_vrs( dfftp%nnr, nspin, doublegrid, kedtau, v%kin_r, vrs )
     !
     ! ... in the US case we have to recompute the self-consistent
     ! ... term in the nonlocal potential
     ! ... PAW: newd contains PAW updates of NL coefficients
     !
     IF (.not. use_gpu) CALL newd()
     IF (      use_gpu) CALL newd_gpu()
     !
  !
END SUBROUTINE qepy_v_of_rho_all
