!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE qepy_electrons_nscf ( printout, exxen)
  !----------------------------------------------------------------------------
  !! This routine is a driver of the self-consistent cycle.
  !! It uses the routine c_bands for computing the bands at fixed
  !! Hamiltonian, the routine sum_band to compute the charge density,
  !! the routine v_of_rho to compute the new potential and the routine
  !! mix_rho to mix input and output charge densities.
  !
  USE kinds,                ONLY : DP
  USE check_stop,           ONLY : check_stop_now, stopped_by_user
  USE io_global,            ONLY : stdout, ionode
  USE cell_base,            ONLY : at, bg, alat, omega, tpiba2
  USE ions_base,            ONLY : zv, nat, nsp, ityp, tau, compute_eextfor, atm
  USE basis,                ONLY : starting_pot
  USE bp,                   ONLY : lelfield
  USE fft_base,             ONLY : dfftp
  USE gvect,                ONLY : ngm, gstart, g, gg, gcutm
  USE gvecs,                ONLY : doublegrid, ngms
  USE klist,                ONLY : xk, wk, nelec, ngk, nks, nkstot, lgauss, &
                                   two_fermi_energies, tot_charge
  USE lsda_mod,             ONLY : lsda, nspin, magtot, absmag, isk
  USE vlocal,               ONLY : strf
  USE wvfct,                ONLY : nbnd, et
  USE gvecw,                ONLY : ecutwfc
  USE ener,                 ONLY : etot, hwf_energy, eband, deband, ehart, &
                                   vtxc, etxc, etxcc, ewld, demet, epaw, &
                                   elondon, edftd3, ef_up, ef_dw, exdm, ef
  USE scf,                  ONLY : scf_type, scf_type_COPY, bcast_scf_type,&
                                   create_scf_type, destroy_scf_type, &
                                   open_mix_file, close_mix_file, &
                                   rho, rho_core, rhog_core, v, vltot, vrs, &
                                   kedtau, vnew
  USE control_flags,        ONLY : mixing_beta, tr2, ethr, niter, nmix, &
                                   iprint, conv_elec, &
                                   restart, io_level, do_makov_payne,  &
                                   gamma_only, iverbosity, textfor,     &
                                   llondon, ldftd3, scf_must_converge, lxdm, ts_vdw
  USE control_flags,        ONLY : n_scf_steps, scf_error

  USE io_files,             ONLY : iunmix, output_drho
  USE ldaU,                 ONLY : eth, Hubbard_U, Hubbard_lmax, &
                                   niter_with_fixed_ns, lda_plus_u
  USE extfield,             ONLY : tefield, etotefield, gate, etotgatefield !TB
  USE noncollin_module,     ONLY : noncolin, magtot_nc, i_cons,  bfield, &
                                   lambda, report
  USE spin_orb,             ONLY : domag
  USE io_rho_xml,           ONLY : write_scf
  USE uspp,                 ONLY : okvan
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp_pools,             ONLY : root_pool, me_pool, my_pool_id, &
                                   inter_pool_comm, intra_pool_comm
  USE mp,                   ONLY : mp_sum, mp_bcast
  !
  USE london_module,        ONLY : energy_london
  USE dftd3_api,            ONLY : dftd3_pbc_dispersion, &
                                   dftd3_init, dftd3_set_functional, &
                                   get_atomic_number, dftd3_input, &
                                   dftd3_calc
  USE dftd3_qe,             ONLY : dftd3, dftd3_in, energy_dftd3
  USE xdm_module,           ONLY : energy_xdm
  USE tsvdw_module,         ONLY : EtsvdW
  !
  USE paw_variables,        ONLY : okpaw, ddd_paw, total_core_energy, only_paw
  USE paw_onecenter,        ONLY : PAW_potential
  USE paw_symmetry,         ONLY : PAW_symmetrize_ddd
  USE dfunct,               ONLY : newd
  USE esm,                  ONLY : do_comp_esm, esm_printpot, esm_ewald
  USE fcp_variables,        ONLY : lfcpopt, lfcpdyn
  USE wrappers,             ONLY : memstat
  !
  USE plugin_variables,     ONLY : plugin_etot
  !
  USE qepy_common,          ONLY : embed
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT (IN) :: printout
  !! * If printout>0, prints on output the total energy;
  !! * if printout>1, also prints decomposition into energy contributions.
  REAL(DP),INTENT (IN) :: exxen
  !! current estimate of the exchange energy
  !
  ! ... local variables
  !
  REAL(DP),save :: dr2
  !! the norm of the diffence between potential
  REAL(DP) :: charge
  !! the total charge
  REAL(DP) :: deband_hwf
  !! deband for the Harris-Weinert-Foulkes functional
  REAL(DP) :: mag
  !! local magnetization
  INTEGER :: i
  !! counter on polarization
  INTEGER :: idum
  !! dummy counter on iterations
  INTEGER,save :: iter
  !! counter on iterations
  INTEGER :: ios, kilobytes
  !
  REAL(DP) :: tr2_min
  !! estimated error on energy coming from diagonalization
  REAL(DP),save :: descf
  !! correction for variational energy
  REAL(DP) :: en_el=0.0_DP
  !! electric field contribution to the total energy
  REAL(DP) :: eext=0.0_DP
  !! external forces contribution to the total energy
  LOGICAL :: first, exst
  !! auxiliary variables for calculating and storing temporary copies of
  !! the charge density and of the HXC-potential
  !
  TYPE(scf_type),save :: rhoin
  !! used to store rho_in of current/next iteration
  !
  ! ... external functions
  !
  REAL(DP), EXTERNAL :: ewald, get_clock
  REAL(DP) :: etot_cmp_paw(nat,2,2)
  ! 
  REAL(DP) :: latvecs(3,3)
  !! auxiliary variables for grimme-d3
  INTEGER:: atnum(1:nat), na
  !! auxiliary variables for grimme-d3
  !
  INTEGER:: its
  !!
  REAL(DP) :: mixing_beta_new
  !qepy <-- add descf
  TYPE(scf_type),save :: rho_prev
  !! save the last step unmix density for descf
  !
  LOGICAL :: add_descf
  !! add the descf to the total energy for last step
  !
  !
  add_descf = .FALSE.
  if ( embed%ldescf .and. embed%iterative .and. (.not. embed%initial) ) add_descf = .TRUE.
  if (embed%finish) goto 10
  !!! If we change some parts to functions will make the code clean.
  !!! But to keep the structure of the code, we add many goto functions.
  !qepy -->

  if (embed%initial) then
  embed%initial = .FALSE.
  iter = 0
  dr2  = 0.0_dp
  descf = 0.0_dp
  IF ( restart ) CALL restart_in_electrons( iter, dr2, ethr, et )
  !end if
  !
  WRITE( stdout, 9000 ) get_clock( 'PWSCF' )
  !
  CALL memstat( kilobytes )
  IF ( kilobytes > 0 ) WRITE( stdout, 9001 ) kilobytes/1000.0
  !
  CALL start_clock( 'electrons' )
  !
  FLUSH( stdout )
  !
  ! ... calculates the ewald contribution to total energy
  !
  !if (embed%initial) then
  if (embed%lewald) then
  IF ( do_comp_esm ) THEN
     ewld = esm_ewald()
  ELSE
     ewld = ewald( alat, nat, nsp, ityp, zv, at, bg, tau, &
                omega, g, gg, ngm, gcutm, gstart, gamma_only, strf )
  ENDIF
  endif
  if (iand(embed%exttype,1) == 1) then
     call qepy_setlocal()
  endif
  !
  IF ( llondon ) THEN
     elondon = energy_london( alat , nat , ityp , at ,bg , tau )
  ELSE
     elondon = 0.d0
  ENDIF
  !
  ! Grimme-D3 correction to the energy
  !
  IF (ldftd3) THEN
     latvecs(:,:)=at(:,:)*alat
     tau(:,:)=tau(:,:)*alat
     DO na = 1, nat
        atnum(na) = get_atomic_number(TRIM(atm(ityp(na))))
     ENDDO
     call dftd3_pbc_dispersion(dftd3,tau,atnum,latvecs,energy_dftd3)
     edftd3=energy_dftd3*2.d0
     tau(:,:)=tau(:,:)/alat
  ELSE
     edftd3= 0.0
  ENDIF
  !
  !
  WRITE( stdout, 9002 )
  FLUSH( stdout )
  !
  else
  dr2 = embed%dnorm
  end if ! if (embed%initial)
  !
  CALL qepy_v_of_rho_all( rho, rho_core, rhog_core, &
     ehart, etxc, vtxc, eth, etotefield, charge, v)

  IF ( check_stop_now() ) THEN
     conv_elec=.FALSE.
     CALL save_in_electrons (iter, dr2, ethr, et )
     GO TO 10
  ENDIF

  iter = iter + 1
  !
  WRITE( stdout, 9010 ) iter, ecutwfc, mixing_beta
  !
  FLUSH( stdout )
  !
  ! ... Convergence threshold for iterative diagonalization is
  ! ... automatically updated during self consistency
  !
  IF ( iter > 1 ) THEN
     !
     IF ( iter == 2 ) ethr = 1.D-2
     ethr = MIN( ethr, 0.1D0*dr2 / MAX( 1.D0, nelec ) )
     ethr = MIN( ethr, embed%diag_conv)
     ! ... do not allow convergence threshold to become too small:
     ! ... iterative diagonalization may become unstable
     ethr = MAX( ethr, 1.D-13 )
     !
  ENDIF
  !
  first = ( iter == 1 )
  !
  if (first) ethr = MAX( ethr, 1.D-6 )
  !
  ! ... tr2_min is set to an estimate of the error on the energy
  ! ... due to diagonalization - used only for the first scf iteration
  !
  tr2_min = 0.D0
  !
  IF ( first ) tr2_min = ethr*MAX( 1.D0, nelec ) 
  !
  ! ... diagonalization of the KS hamiltonian
  !
  IF ( lelfield ) THEN
     CALL c_bands_efield( iter )
  ELSE
     CALL c_bands( iter )
  ENDIF
  !
  IF ( stopped_by_user ) THEN
     conv_elec=.FALSE.
     CALL save_in_electrons( iter-1, dr2, ethr, et )
     GO TO 10
  ENDIF
  !
  ! ... xk, wk, isk, et, wg are distributed across pools;
  ! ... the first node has a complete copy of xk, wk, isk,
  ! ... while eigenvalues et and weights wg must be
  ! ... explicitly collected to the first node
  ! ... this is done here for et, in sum_band for wg
  !
  CALL poolrecover( et, nbnd, nkstot, nks )
  !
  ! ... the new density is computed here. For PAW:
  ! ... sum_band computes new becsum (stored in uspp modules)
  ! ... and a subtly different copy in rho%bec (scf module)
  !
  !CALL sum_band()
  RETURN

10  FLUSH( stdout )
  !qepy <--
  n_scf_steps = iter
  scf_error = dr2
  !qepy -->
  !
  ! ... exiting: write (unless disabled) the charge density to file
  ! ... (also write ldaU ns coefficients and PAW becsum)
  !
  IF ( io_level > -1 ) CALL write_scf( rho, nspin )
  !
  ! ... delete mixing info if converged, keep it if not
  !
  IF ( embed%finish) conv_elec = .true.
  IF ( output_drho /= ' ' ) CALL remove_atomic_rho()
  CALL stop_clock( 'electrons' )
  !
  !qepy <-- reset embed
  embed%initial = .TRUE.
  embed%finish = .FALSE.
  embed%mix_coef = -1.0
  !qepy -->
  IF ( conv_elec ) THEN
     WRITE( stdout, 9101 )
     WRITE( stdout, 9110 ) iter
  ELSE
     WRITE( stdout, 9120 ) iter
  ENDIF
  RETURN
  !
  ! ... formats
  !
9000 FORMAT(/'     total cpu time spent up to now is ',F10.1,' secs' )
9001 FORMAT(/'     per-process dynamical memory: ',f7.1,' Mb' )
9002 FORMAT(/'     Self-consistent Calculation' )
9010 FORMAT(/'     iteration #',I3,'     ecut=', F9.2,' Ry',5X,'beta=',F5.2 )
9050 FORMAT(/'     WARNING: integrated charge=',F15.8,', expected=',F15.8 )
9101 FORMAT(/'     End of self-consistent calculation' )
9110 FORMAT(/'     convergence has been achieved in ',i3,' iterations' )
9120 FORMAT(/'     convergence NOT achieved after ',i3,' iterations: stopping' )
  !
END SUBROUTINE qepy_electrons_nscf
