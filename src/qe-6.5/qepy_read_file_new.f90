!
! Copyright (C) 2016-2019 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE qepy_read_file()
  !----------------------------------------------------------------------------
  !
  ! Wrapper routine, for backwards compatibility
  !
  USE io_global,        ONLY : stdout
  USE control_flags,    ONLY : io_level
  USE buffers,          ONLY : open_buffer, close_buffer, save_buffer
  USE io_files,         ONLY : nwordwfc, iunwfc, restart_dir
  USE wvfct,            ONLY : nbnd, npwx
  USE noncollin_module, ONLY : npol
  USE klist,            ONLY : nks
  USE wavefunctions,    ONLY : evc
  USE pw_restart_new,   ONLY : read_collected_wfc
  !
  IMPLICIT NONE
  !
  INTEGER :: ik
  LOGICAL :: exst, wfc_is_collected
  !
  wfc_is_collected = .true.
  CALL qepy_read_file_new( wfc_is_collected )
  !
  ! ... Open unit iunwfc, for Kohn-Sham orbitals - we assume that wfcs
  ! ... have been written to tmp_dir, not to a different directory!
  ! ... io_level = 1 so that a real file is opened
  !
  nwordwfc = nbnd*npwx*npol
  io_level = 1
  CALL open_buffer ( iunwfc, 'wfc', nwordwfc, io_level, exst )
  !
  ! ... read wavefunctions in collected format, write them to file
  !
  IF ( wfc_is_collected ) THEN
     !
     WRITE( stdout, '(5x,A)') &
          'Reading collected, re-writing distributed wavefunctions'
     DO ik = 1, nks
        CALL read_collected_wfc ( restart_dir(), ik, evc )
        CALL save_buffer ( evc, nwordwfc, iunwfc, ik )
     END DO
     !
  ELSE
     WRITE( stdout, '(5x,A)') &
          'read_file: Wavefunctions in collected format not available'
  END IF
  !
  !CALL close_buffer  ( iunwfc, 'KEEP' )
  !
END SUBROUTINE qepy_read_file
!
!----------------------------------------------------------------------------
SUBROUTINE qepy_read_file_new ( needwf )
  !----------------------------------------------------------------------------
  !
  ! Reads xml data file produced by pw.x or cp.x, performs initializations
  ! related to the contents of the xml file
  ! If needwf=.t. performs wavefunction-related initialization as well
  ! Does not read wfcs but returns in "wfc_is_collected" info on the wfc file
  !
  USE io_global,      ONLY : stdout
  USE io_files,       ONLY : nwordwfc, iunwfc, wfc_dir, tmp_dir, restart_dir
  USE gvect,          ONLY : ngm, g
  USE gvecw,          ONLY : gcutw
  USE klist,          ONLY : nkstot, nks, xk, wk
  USE lsda_mod,       ONLY : isk
  USE wvfct,          ONLY : nbnd, et, wg
  USE qepy_pw_restart_new, ONLY : qepy_read_xml_file
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(INOUT) :: needwf
  !
  LOGICAL :: wfc_is_collected
  !
  WRITE( stdout, '(/,5x,A)') &
       'Reading xml data from directory:', TRIM( restart_dir() )
  !
  ! ... Read the contents of the xml data file
  !
  CALL qepy_read_xml_file ( wfc_is_collected )
  !
  ! ... more initializations: pseudopotentials / G-vectors / FFT arrays /
  ! ... charge density / potential / ... , but not KS orbitals
  !
  CALL qepy_post_xml_init ( )
  !
  IF ( needwf ) THEN
     IF ( .NOT. wfc_is_collected ) WRITE( stdout, '(5x,A)') &
          'qepy_read_file_new: Wavefunctions not in collected format?!?'
     !
     ! ... initialization of KS orbitals
     !
     wfc_dir = tmp_dir ! this is likely obsolete and no longer used
     !
     ! ... distribute across pools k-points and related variables.
     ! ... nks is defined by the following routine as the number 
     ! ... of k-points in the current pool
     !
     CALL divide_et_impera( nkstot, xk, wk, isk, nks )
     CALL poolscatter( nbnd, nkstot, et, nks, et )
     CALL poolscatter( nbnd, nkstot, wg, nks, wg )
     !
     ! ... allocate_wfc_k also computes no. of plane waves and k+G indices
     ! ... FIXME: the latter should be read from file, not recomputed
     !
     CALL allocate_wfc_k()
     !
  END IF
  needwf = wfc_is_collected
  !
END SUBROUTINE qepy_read_file_new
!----------------------------------------------------------------------------
SUBROUTINE qepy_post_xml_init (  )
  !----------------------------------------------------------------------------
  !
  ! ... Various initializations needed to start a calculation:
  ! ... pseudopotentials, G vectors, FFT arrays, rho, potential
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE uspp_param,           ONLY : upf
  USE read_pseudo_mod,      ONLY : readpp
  USE uspp,                 ONLY : becsum
  USE paw_variables,        ONLY : okpaw, ddd_PAW
  USE paw_init,             ONLY : paw_init_onecenter, allocate_paw_internals
  USE paw_onecenter,        ONLY : paw_potential
  USE dfunct,               ONLY : newd
  USE funct,                ONLY : get_inlc, get_dft_name
  USE ldaU,                 ONLY : lda_plus_u, eth, init_lda_plus_u, U_projection
  USE esm,                  ONLY : do_comp_esm, esm_init
  USE Coul_cut_2D,          ONLY : do_cutoff_2D, cutoff_fact 
  USE ions_base,            ONLY : nat, nsp, tau, ityp
  USE recvec_subs,          ONLY : ggen, ggens
  USE gvect,                ONLY : gg, ngm, g, gcutm, mill, ngm_g, ig_l2g, &
                                   eigts1, eigts2, eigts3, gstart, gshells
  USE gvecs,                ONLY : ngms, gcutms 
  USE fft_rho,              ONLY : rho_g2r
  USE fft_base,             ONLY : dfftp, dffts
  USE scf,                  ONLY : rho, rho_core, rhog_core, v
  USE io_rho_xml,           ONLY : read_scf
  USE vlocal,               ONLY : strf
  USE control_flags,        ONLY : gamma_only
  USE control_flags,        ONLY : ts_vdw, tqr, tq_smoothing, tbeta_smoothing
  USE cellmd,               ONLY : cell_factor, lmovecell
  USE wvfct,                ONLY : nbnd, nbndx, et, wg
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module,     ONLY : noncolin
  USE spin_orb,             ONLY : lspinorb
  USE cell_base,            ONLY : at, bg, set_h_ainv
  USE symm_base,            ONLY : d1, d2, d3
  USE realus,               ONLY : betapointlist, generate_qpointlist, &
                                   init_realspace_vars,real_space
  !
  !qepy --> import
  USE control_flags,        ONLY : mixing_beta, tr2, ethr, niter, nmix, &
                                   iprint, conv_elec, &
                                   restart, io_level, do_makov_payne,  &
                                   iverbosity, textfor,     &
                                   llondon, ldftd3, scf_must_converge, lxdm
  
  USE london_module,        ONLY : energy_london, init_london, C6_ij
  USE dftd3_api,            ONLY : dftd3_pbc_dispersion, &
                                   dftd3_init, dftd3_set_functional, &
                                   get_atomic_number, dftd3_input, &
                                   dftd3_calc
  USE dftd3_qe,             ONLY : dftd3, dftd3_in, energy_dftd3, dftd3_xc, dftd3_printout
  USE xdm_module,           ONLY : energy_xdm, init_xdm
  USE input_parameters,     ONLY : dftd3_threebody, dftd3_version
  USE funct,                ONLY : get_dft_short
  USE tsvdw_module,         ONLY : tsvdw_initialize
  !qepy <-- import
  IMPLICIT NONE
  !
  INTEGER  :: inlc
  REAL(DP) :: ehart, etxc, vtxc, etotefield, charge
  CHARACTER(LEN=20) :: dft_name
  !qepy --> variablesi
  CHARACTER(LEN=256):: dft_
  REAL (DP), EXTERNAL :: get_clock
  !qepy <-- variables
  !
  ! ... set G cutoffs and cell factor (FIXME: from setup.f90?)
  !
  CALL set_gcut()
  if (cell_factor == 0.d0) cell_factor = 1.D0
  nbndx = nbnd
  !
  ! ... read pseudopotentials
  ! ... the following call prevents readpp from setting dft from PP files
  !
  dft_name = get_dft_name ()
  CALL readpp ( dft_name )
  !
  ! ... misc PP initialization (from setup.f90)
  !
  okpaw = ANY ( upf(1:nsp)%tpawp )
  IF ( .NOT. lspinorb ) CALL average_pp ( nsp )
  !! average_pp must be called before init_lda_plus_u
  IF ( lda_plus_u ) CALL init_lda_plus_u ( upf(1:nsp)%psd, noncolin )
  !
  ! ... allocate memory for G- and R-space fft arrays (from init_run.f90)
  !
  CALL pre_init()
  ! NB: data_structure uses k-points to compute gkcut
  CALL data_structure ( gamma_only )
  CALL allocate_fft()
  CALL ggen ( dfftp, gamma_only, at, bg, gcutm, ngm_g, ngm, &
       g, gg, mill, ig_l2g, gstart ) 
  CALL ggens( dffts, gamma_only, at, g, gg, mill, gcutms, ngms ) 
  CALL gshells ( lmovecell ) 
  !
  IF (do_comp_esm) CALL esm_init()
  IF (do_cutoff_2D) CALL cutoff_fact()
  !
  ! ... allocate the potentials
  !
  CALL allocate_locpot()
  CALL allocate_nlpot()
  IF (okpaw) THEN
     CALL allocate_paw_internals()
     CALL paw_init_onecenter()
     CALL d_matrix(d1,d2,d3)
  ENDIF
  !
  ! ... read the charge density in G-space
  !
  CALL read_scf( rho, nspin, gamma_only )
  !
  ! ... bring the charge density to real space
  !
  CALL rho_g2r ( dfftp, rho%of_g, rho%of_r )
  !
  ! ... re-compute the local part of the pseudopotential vltot and
  ! ... the core correction charge (if any) - from hinit0.f90
  !
  CALL init_vloc()
  IF (tbeta_smoothing) CALL init_us_b0()
  IF (tq_smoothing) CALL init_us_0()
  CALL init_us_1()
  IF ( lda_plus_U .AND. ( U_projection == 'pseudo' ) ) CALL init_q_aeps()
  CALL init_at_1()
  !
  CALL struc_fact( nat, tau, nsp, ityp, ngm, g, bg, dfftp%nr1, dfftp%nr2,&
                   dfftp%nr3, strf, eigts1, eigts2, eigts3 )
  CALL setlocal()
  CALL set_rhoc()
  !
  ! ... for real-space PP's
  !
  IF ( tqr ) CALL generate_qpointlist()
  IF (real_space ) THEN
     CALL betapointlist()
     CALL init_realspace_vars()
     WRITE (stdout,'(5X,"Real space initialisation completed")')    
  ENDIF
  !
  ! ... recalculate the potential - FIXME: couldn't make ts-vdw work
  !
  !qepy --> init vdw
  IF ( ts_vdw) THEN
      CALL tsvdw_initialize()
      CALL set_h_ainv()
     !CALL infomsg('qepy_read_file_new','*** vdW-TS term will be missing in potential ***')
     !ts_vdw = .false.
  END IF
  !
  !IF ( llondon ) THEN
  !   IF ( .NOT. ALLOCATED(C6_ij)) CALL init_london()
  !ENDIF
  !IF (lxdm) THEN
  !   IF ( get_clock('init_xdm')< 0.0_DP) CALL init_xdm()
  !ENDIF
  !IF (ldftd3) !   IF ( .NOT. ALLOCATED(dftd3%r0ab)) THEN
  IF (llondon) CALL init_london()
  IF (lxdm) CALL init_xdm()
  IF (ldftd3) THEN
     if (dftd3_version==2) dftd3_threebody=.false.
     dftd3_in%threebody = dftd3_threebody
     CALL dftd3_init(dftd3, dftd3_in)
     CALL dftd3_printout(dftd3, dftd3_in)
     dft_ = get_dft_short( )
     dft_ = dftd3_xc ( dft_ )
     CALL dftd3_set_functional(dftd3, func=dft_, version=dftd3_version,tz=.false.)
  ENDIF
  !qepy <-- init
  CALL v_of_rho( rho, rho_core, rhog_core, &
       ehart, etxc, vtxc, eth, etotefield, charge, v )
  !
  ! ... More PAW and USPP initializations
  !
  IF (okpaw) THEN
     becsum = rho%bec
     CALL PAW_potential(rho%bec, ddd_PAW)
  ENDIF 
  CALL newd()
  !
  RETURN
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE set_gcut()
      !------------------------------------------------------------------------
      !
      USE constants, ONLY : pi, eps8
      USE cell_base, ONLY : alat, tpiba, tpiba2
      USE gvect,     ONLY : ecutrho, gcutm
      USE gvecs,     ONLY : gcutms, dual, doublegrid
      USE gvecw,     ONLY : gcutw, ecutwfc
      !
      !
      ! ... Set the units in real and reciprocal space
      !
      tpiba  = 2.D0 * pi / alat
      tpiba2 = tpiba**2
      !
      ! ... Compute the cut-off of the G vectors
      !
      gcutw =        ecutwfc / tpiba2
      gcutm = dual * ecutwfc / tpiba2
      ecutrho=dual * ecutwfc
      !
      doublegrid = ( dual > 4.0_dp + eps8 )
      IF ( doublegrid ) THEN
         gcutms = 4.D0 * ecutwfc / tpiba2
      ELSE
         gcutms = gcutm
      END IF
      !
    END SUBROUTINE set_gcut
    !
  END SUBROUTINE qepy_post_xml_init
