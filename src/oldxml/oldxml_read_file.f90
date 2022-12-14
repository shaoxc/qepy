!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#if defined(__OLDXML)
!----------------------------------------------------------------------------
! TB
! included allocation of the force field of the gate, search for 'TB'
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
SUBROUTINE oldxml_read_file()
  !----------------------------------------------------------------------------
  !
  ! Wrapper routine, for compatibility
  !
  USE io_files,             ONLY : nwordwfc, iunwfc, prefix, tmp_dir, wfc_dir, diropn
  USE io_global,            ONLY : stdout, ionode
  USE buffers,              ONLY : open_buffer, close_buffer
  USE wvfct,                ONLY : nbnd, npwx
  USE noncollin_module,     ONLY : npol
  USE paw_variables,        ONLY : okpaw, ddd_PAW
  USE paw_onecenter,        ONLY : paw_potential
  USE uspp,                 ONLY : becsum
  USE scf,                  ONLY : rho
  USE realus,               ONLY : betapointlist, &
                                   init_realspace_vars,real_space
  USE dfunct,               ONLY : newd
  USE ldaU,                 ONLY : lda_plus_u, U_projection
  USE oldxml_pw_restart,    ONLY : pw_readfile
  USE control_flags,        ONLY : io_level
  USE klist,                ONLY : init_igk
  USE gvect,                ONLY : ngm, g
  USE gvecw,                ONLY : gcutw
#if defined (__HDF5)
  USE hdf5_qe
#endif
  USE klist,                ONLY : nks
  USE wavefunctions,        ONLY : evc
  USE buffers,              ONLY : get_buffer, save_buffer
  !
  IMPLICIT NONE 
  INTEGER :: ierr
  LOGICAL :: exst, exst_mem, exst_file, opnd_file
  CHARACTER( 256 )  :: dirname
  INTEGER  :: ik
  INTEGER  :: file_size
  !
  !
  ierr = 0 
  !
  ! ... Read the contents of the xml data file
  !
  IF ( ionode ) WRITE( stdout, '(/,5x,A,/,5x,A)') &
     'Reading data from directory:', TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
#if defined(__HDF5)
  CALL initialize_hdf5()
#endif
  !
  CALL oldxml_read_xml_file ( )
  !
  ! ... Open unit iunwfc, for Kohn-Sham orbitals - we assume that wfcs
  ! ... have been written to tmp_dir, not to a different directory!
  ! ... io_level = 1 so that a real file is opened
  !
  wfc_dir = tmp_dir
  nwordwfc = nbnd*npwx*npol
  io_level = 1
  !CALL open_buffer ( iunwfc, 'wfc', nwordwfc, io_level, exst )
  CALL open_buffer( iunwfc, 'wfc', nwordwfc, io_level, exst_mem, exst_file )
  !
  ! ... Allocate and compute k+G indices and number of plane waves
  ! ... FIXME: should be read from file, not re-computed
  !
  CALL init_igk ( npwx, ngm, g, gcutw ) 
  !
  !CALL pw_readfile( 'wave', ierr )
  !IF ( nks == 1 ) CALL get_buffer( evc, nwordwfc, iunwfc, 1 )
  !-----------------------------------------------------------------------
  !CALL open_buffer( iunwfc, 'wfc', nwordwfc, io_level, exst_mem, exst_file )
  !
  ! ... now the various possible wavefunction initializations
  ! ... first a check: is "tmp_dir"/"prefix".wfc found on disk?
  !
  IF ( .NOT. exst_file) THEN
     !
     ! ... "tmp_dir"/"prefix".wfc not found on disk: try to read
     ! ... wavefunctions in "collected" format from "prefix".save/, 
     ! ... rewrite them (in pw_readfile) using the internal format
     !
     ierr = 1
     CALL pw_readfile( 'wave', ierr )
     IF ( ierr > 0 ) THEN
        WRITE( stdout, '(5x,A)') &
             'oldxml_read_file: Wavefunctions in collected format not available'
        !CALL errore( 'read_file ', 'file not found', ierr )
     !END IF
     ELSE
     !
     ! ... workaround: with k-point parallelization and 1 k-point per pool,
     ! ... pw_readfile does not leave evc properly initialized on all pools
     !
     inquire (unit = iunwfc, size = file_size)
     !qepy --> skip empty file
     if (file_size /= 0) then
     IF ( nks == 1 ) CALL get_buffer( evc, nwordwfc, iunwfc, 1 )
     else
        WRITE( stdout, '(5x,A)') 'oldxml_read_file: Wavefunctions is empty!'
     endif
     !qepy <-- skip empty file
     END IF
     !
  ELSE
     !
     ! ... wavefunctions are read from file (or buffer) in routine 
     ! ... c_bands, but not if there is a single k-point. In such
     ! ... a case, we read wavefunctions (directly from file in 
     ! ... order to avoid a useless buffer allocation) here
     !
     IF ( nks == 1 ) THEN
         inquire (unit = iunwfc, opened = opnd_file)
         if (.not.opnd_file) CALL diropn( iunwfc, 'wfc', 2*nwordwfc, exst )
         inquire (unit = iunwfc, size = file_size)
         !qepy --> skip empty file
         if (file_size /= 0) then
         CALL davcio ( evc, 2*nwordwfc, iunwfc, nks, -1 )
         else
            WRITE( stdout, '(5x,A)') 'oldxml_read_file: Wavefunctions is empty!'
         endif
         !qepy <-- skip empty file
         if(.not.opnd_file) CLOSE ( UNIT=iunwfc, STATUS='keep' )
     END IF
     !
  END IF
  !-----------------------------------------------------------------------
  !
  ! ... Assorted initialization: pseudopotentials, PAW
  ! ... Not sure which ones (if any) should be done here
  !
  CALL init_us_1()
  !
  IF (lda_plus_u .AND. (U_projection == 'pseudo')) CALL init_q_aeps()
  !
  IF (okpaw) THEN
     becsum = rho%bec
     CALL PAW_potential(rho%bec, ddd_PAW)
  ENDIF 
  !
  IF ( real_space ) THEN
    CALL betapointlist()
    CALL init_realspace_vars()
    IF( ionode ) WRITE(stdout,'(5x,"Real space initialisation completed")')
  ENDIF
  CALL newd()
  !-----------------------------------------------------------------------
  !IF ( wfc_is_collected ) THEN
  !   !
  !   WRITE( stdout, '(5x,A)') &
  !        'Reading collected, re-writing distributed wavefunctions'
  !   DO ik = 1, nks
  !      CALL read_collected_wfc ( restart_dir(), ik, evc )
  !      CALL save_buffer ( evc, nwordwfc, iunwfc, ik )
  !   END DO
  !   !
  !ELSE
  !   WRITE( stdout, '(5x,A)') &
  !        'read_file: Wavefunctions in collected format not available'
  !END IF
  !-----------------------------------------------------------------------
  !
  !CALL close_buffer  ( iunwfc, 'KEEP' )
  !
END SUBROUTINE oldxml_read_file
!
SUBROUTINE oldxml_read_xml_file()
  ! wrapper routine to call the default behavior
  call oldxml_read_xml_file_internal(.true.)
END SUBROUTINE oldxml_read_xml_file

SUBROUTINE oldxml_read_xml_file_nobs()
  ! wrapper routine to load everything except for the band structure
  call oldxml_read_xml_file_internal(.false.)
END SUBROUTINE oldxml_read_xml_file_nobs

!----------------------------------------------------------------------------
SUBROUTINE oldxml_read_xml_file_internal(withbs)
  !----------------------------------------------------------------------------
  !
  ! ... This routine allocates space for all quantities already computed
  ! ... in the pwscf program and reads them from the data file.
  ! ... All quantities that are initialized in subroutine "setup" when
  ! ... starting from scratch should be initialized here when restarting
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, nsp, ityp, tau, if_pos, extfor
  USE cell_base,            ONLY : tpiba2, alat,omega, at, bg, ibrav
  USE force_mod,            ONLY : force
  USE klist,                ONLY : nkstot, nks, xk, wk
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE wvfct,                ONLY : nbnd, nbndx, et, wg
  USE symm_base,            ONLY : irt, d1, d2, d3, checkallsym, nsym
  USE extfield,             ONLY : forcefield, tefield, gate, forcegate
  USE cellmd,               ONLY : cell_factor, lmovecell
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE fft_types,            ONLY : fft_type_allocate
  USE recvec_subs,          ONLY : ggen, ggens
  USE gvect,                ONLY : gg, ngm, g, gcutm, mill, ngm_g, ig_l2g, &
                                   eigts1, eigts2, eigts3, gstart, gshells
  USE Coul_cut_2D,          ONLY : do_cutoff_2D, cutoff_fact
  USE fft_base,             ONLY : dfftp, dffts
  USE gvecs,                ONLY : ngms, gcutms 
  USE spin_orb,             ONLY : lspinorb, domag
  USE scf,                  ONLY : rho, rho_core, rhog_core, v
  USE wavefunctions,        ONLY : psic
  USE vlocal,               ONLY : strf
  USE io_files,             ONLY : tmp_dir, prefix, iunpun, nwordwfc, iunwfc
  USE noncollin_module,     ONLY : noncolin, npol, nspin_lsda, nspin_mag, nspin_gga
  USE oldxml_pw_restart,           ONLY : pw_readfile, pp_check_file
  USE oldxml_io_rho_xml,           ONLY : read_scf
  USE read_pseudo_mod,      ONLY : readpp
  USE uspp,                 ONLY : becsum
  USE uspp_param,           ONLY : upf
  USE paw_variables,        ONLY : okpaw, ddd_PAW
  USE paw_init,             ONLY : paw_init_onecenter, allocate_paw_internals
  USE ldaU,                 ONLY : lda_plus_u, eth, init_lda_plus_u
  USE control_flags,        ONLY : gamma_only
  USE funct,                ONLY : get_inlc, get_dft_name
  !USE kernel_table,         ONLY : initialize_kernel_table
  USE esm,                  ONLY : do_comp_esm, esm_init
  USE mp_bands,             ONLY : intra_bgrp_comm, nyfft
  USE spin_orb,             ONLY : lforcet
  USE io_files,             ONLY : postfix
  USE oldxml_xml_io_base,          ONLY : read_rho
  USE wvfct,                ONLY : nbnd, npwx
  !
  IMPLICIT NONE

  ! Used to specify whether to read the band structure (files 
  ! K??????/eigenval.xml), so one can skip it if not needed by
  ! the post-processing tool. 
  ! Set to True for the 'default' behavior of reading these files.
  LOGICAL :: withbs

  INTEGER  :: i, is, ik, ibnd, nb, nt, ios, isym, ierr, inlc
  REAL(DP) :: rdum(1,1), ehart, etxc, vtxc, etotefield, charge
  REAL(DP) :: sr(3,3,48)
  CHARACTER(LEN=20) dft_name
  CHARACTER(LEN=256)    :: dirname, filename
  !
  !
  ! ... first we get the version of the qexml file
  !     if not already read
  CALL pw_readfile( 'header', ierr )
  CALL errore( 'read_xml_file ', 'unable to determine qexml version', ABS(ierr) )
  !
  ! ... then we check if the file can be used for post-processing
  !
  IF ( .NOT. pp_check_file() ) CALL infomsg( 'read_xml_file', &
               & 'file ' // TRIM( tmp_dir ) // TRIM( prefix ) &
               & // '.save not guaranteed to be safe for post-processing' )
  !
  ! ... here we read the variables that dimension the system
  ! ... in parallel execution, only root proc reads the file
  ! ... and then broadcasts the values to all other procs
  !
  CALL pw_readfile( 'reset', ierr )
  CALL pw_readfile( 'dim',   ierr )
  CALL errore( 'read_xml_file ', 'problem reading file ' // &
             & TRIM( tmp_dir ) // TRIM( prefix ) // '.save', ierr )
  !
  ! ... allocate space for atomic positions, symmetries, forces
  !
  IF ( nat < 0 ) CALL errore( 'read_xml_file', 'wrong number of atoms', 1 )
  !
  ! ... allocation
  !
  ALLOCATE( ityp( nat ) )
  ALLOCATE( tau(    3, nat ) )
  ALLOCATE( if_pos( 3, nat ) )
  ALLOCATE( force(  3, nat ) )
  ALLOCATE( extfor(  3, nat ) )
  !
  IF ( tefield ) ALLOCATE( forcefield( 3, nat ) )
  IF ( gate ) ALLOCATE( forcegate( 3, nat ) ) ! TB
  !
  ALLOCATE( irt( 48, nat ) )
  !
  CALL oldxml_set_dimensions()
  CALL fft_type_allocate ( dfftp, at, bg, gcutm, intra_bgrp_comm, nyfft=nyfft )
  CALL fft_type_allocate ( dffts, at, bg, gcutms, intra_bgrp_comm, nyfft=nyfft )
  !
  ! ... check whether LSDA
  !
  IF ( lsda ) THEN
     !
     nspin = 2
     npol  = 1
     !
  ELSE IF ( noncolin ) THEN
     !
     nspin        = 4
     npol         = 2
     current_spin = 1
     !
  ELSE
     !
     nspin        = 1
     npol         = 1
     current_spin = 1
     !
  END IF
  !
  if (cell_factor == 0.d0) cell_factor = 1.D0
  !
  ! ... allocate memory for eigenvalues and weights (read from file)
  !
  nbndx = nbnd
  ALLOCATE( et( nbnd, nkstot ) , wg( nbnd, nkstot ) )
  !
  ! ... here we read all the variables defining the system
  !
  IF  ( withbs .EQV. .TRUE. ) THEN  
     CALL pw_readfile( 'nowave', ierr )
  ELSE
     CALL pw_readfile( 'nowavenobs', ierr )
  END IF
  !
  ! ... distribute across pools k-points and related variables.
  ! ... nks is defined by the following routine as the number 
  ! ... of k-points in the current pool
  !
  CALL divide_et_impera( nkstot, xk, wk, isk, nks )
  !
  CALL poolscatter( nbnd, nkstot, et, nks, et )
  CALL poolscatter( nbnd, nkstot, wg, nks, wg )
  !
  ! ... check on symmetry
  !
  IF (nat > 0) CALL checkallsym( nat, tau, ityp )
  !
  !  Set the different spin indices
  !
  nspin_mag  = nspin
  nspin_lsda = nspin
  nspin_gga  = nspin
  IF (nspin==4) THEN
     nspin_lsda=1
     IF (domag) THEN
        nspin_gga=2
     ELSE
        nspin_gga=1
        nspin_mag=1
     ENDIF
  ENDIF
  !
  ! ... read pseudopotentials
  !
  CALL pw_readfile( 'pseudo', ierr )

  dft_name = get_dft_name () ! already set, should not be set again
  CALL readpp ( dft_name )
  !
  ! ... read the vdw kernel table if needed
  !
  inlc = get_inlc()
  !if (inlc > 0 ) then
  !    call initialize_kernel_table(inlc)
  !endif
  !
  okpaw = ANY ( upf(1:nsp)%tpawp )
  !
  IF ( .NOT. lspinorb ) CALL average_pp ( nsp )
  !
  ! ... allocate memory for G- and R-space fft arrays
  !
  CALL pre_init()
  CALL data_structure ( gamma_only )
  CALL allocate_fft()
  CALL ggen ( dfftp, gamma_only, at, bg, gcutm, ngm_g, ngm, &
       g, gg, mill, ig_l2g, gstart ) 
  CALL ggens( dffts, gamma_only, at, g, gg, mill, gcutms, ngms ) 
  IF (do_comp_esm) THEN
    CALL pw_readfile( 'esm', ierr )
    CALL esm_init()
  END IF
  CALL gshells ( lmovecell ) 
  !
  ! ... allocate the potential and wavefunctions
  !
  CALL allocate_locpot()
  CALL allocate_nlpot()
  IF (okpaw) THEN
     CALL allocate_paw_internals()
     CALL paw_init_onecenter()
     CALL d_matrix(d1,d2,d3)
  ENDIF
  !
  IF ( lda_plus_u ) THEN
     CALL init_lda_plus_u ( upf(1:nsp)%psd, noncolin )
  ENDIF
  !
  !CALL allocate_wfc()
  CALL allocate_wfc_k()
  !
  ! ... read the charge density
  !
  CALL read_scf( rho, nspin )
!  IF ( .NOT.lforcet ) THEN
!     CALL read_scf ( rho, nspin, gamma_only )
!  ELSE
!     !
!     ! ... 'force theorem' calculation of MAE: read rho only from previous
!     ! ... lsda calculation, set noncolinear magnetization from angles
!     !
!     dirname = TRIM(tmp_dir) // TRIM (prefix) // postfix
!#if defined __HDF5
!     filename = TRIM(dirname) // 'charge-density.hdf5'
!#else 
!     filename = TRIM(dirname) // 'charge-density.dat'
!#endif
!     CALL read_rho ( dirname, rho%of_r, 2 )
!     CALL nc_magnetization_from_lsda ( dfftp%nnr, nspin, rho%of_r )
!  END IF
  !
  ! ... re-calculate the local part of the pseudopotential vltot
  ! ... and the core correction charge (if any) - This is done here
  ! ... for compatibility with the previous version of read_file
  !
  ! 2D calculations: re-initialize cutoff fact before calculating potentials
  IF(do_cutoff_2D) CALL cutoff_fact()
  !
  CALL init_vloc()
  CALL struc_fact( nat, tau, nsp, ityp, ngm, g, bg, dfftp%nr1, dfftp%nr2, &
                   dfftp%nr3, strf, eigts1, eigts2, eigts3 )
  CALL setlocal()
  CALL set_rhoc()
  !
  ! ... bring rho to G-space
  !
  DO is = 1, nspin
     !
     psic(:) = rho%of_r(:,is)
     CALL fwfft ('Rho', psic, dfftp)
     rho%of_g(:,is) = psic(dfftp%nl(:))
     !
  END DO
  !
  ! ... read info needed for hybrid functionals
  !
  CALL pw_readfile('exx', ierr)
  !
  ! ... recalculate the potential
  !
  CALL v_of_rho( rho, rho_core, rhog_core, &
                 ehart, etxc, vtxc, eth, etotefield, charge, v )
  !
  !
  RETURN
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE oldxml_set_dimensions()
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
    END SUBROUTINE oldxml_set_dimensions
    !
  END SUBROUTINE oldxml_read_xml_file_internal
#endif
