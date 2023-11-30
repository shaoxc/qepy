!
! Copyright (C) 2001-2014 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
SUBROUTINE qepy_tddft_readin(infile)
  !-----------------------------------------------------------------------
  !
  ! ... Read in the tddft input file. The input file consists of a
  ! ... single namelist &inputtddft. See doc/user-manual.pdf for the
  ! ... list of input keywords.
  !
  USE tddft_module
  USE io_files,         ONLY : prefix, tmp_dir  
  USE io_global,        ONLY : ionode
  USE constants,        ONLY : bohr_radius_angs, au_sec
  USE mp_images,        ONLY : my_image_id

  ! -- local variables ---------------------------------------------------
  implicit none
  integer :: ios
  character(len=256), external :: trimcheck
  character(len=80) :: verbosity
  !
  CHARACTER(len=*), INTENT(IN), OPTIONAL  :: infile
  INTEGER, EXTERNAL :: find_free_unit
  INTEGER :: iun, i
  CHARACTER(len=256) :: fstr
  CHARACTER(len=1), EXTERNAL :: lowercase
  !
  namelist /inputtddft/ job, prefix, tmp_dir, conv_threshold, verbosity, &
                        dt, e_strength, e_direction, nstep, nupdate_Dnm, &
                        l_circular_dichroism, l_tddft_restart, max_seconds, &
                        molecule, ehrenfest, isave_rho, wavepacket, wp_pos, &
                        wp_d, wp_ekin

  if (.not. ionode .or. my_image_id > 0) goto 400

  ! define input defult values
  call get_environment_variable( 'ESPRESSO_TMPDIR', tmp_dir ) 
  if (trim(tmp_dir) == ' ') tmp_dir = './scratch/'
  tmp_dir = trimcheck(tmp_dir)
  job          = ''
  prefix       = 'pwscf'
  tmp_dir      = './scratch/'    
  verbosity    = 'low'
  dt           = 2.d0                      ! time step (default: 2 attosecond)
  e_strength   = 0.01d0                    ! impulse electric field strength (default: 0.01/Ang)
  e_direction  = 1                         ! impulse electric field direction: 1-x 2-y 3-z
  conv_threshold = 1.0d-12                 ! convergence threshold    
  nstep        = 1000                      ! total time steps
  nupdate_Dnm  = 1                         ! update USPP Dnm every step
  l_circular_dichroism = .false.
  l_tddft_restart      = .false.
  max_seconds  =  1.d7
  molecule     = .true.

  ehrenfest    = .false.                   ! Ehrenfest dynamics
  isave_rho    = 0                         ! save density in XSF format

  wavepacket   = .false.
  wp_pos(1:3) = 0.d0
  wp_d(1:3) = 0.d0
  wp_ekin = 0.d0
 
  ! read input    
  IF ( PRESENT(infile)) THEN
     iun = find_free_unit()
     OPEN ( UNIT = iun, FILE = infile, FORM = 'FORMATTED', &
        STATUS = 'OLD', IOSTAT = ios)
     IF (ios /= 0) GOTO 200
     DO 
        READ(iun, '(A200)', IOSTAT=ios) fstr
        IF (ios /= 0) GOTO 200
        fstr=ADJUSTL(fstr)
        IF (fstr(1:1)=='&') THEN
           DO i=2, LEN_TRIM(fstr)
              fstr(i:i)=lowercase(fstr(i:i))
           ENDDO
           IF (fstr=='&inputtddft') THEN
              BACKSPACE(iun)
              EXIT
           ENDIF
        ENDIF
     ENDDO
     READ( iun, inputtddft, ERR = 200, IOSTAT = ios )
     CLOSE( iun )
  ELSE
  call input_from_file()
  read( 5, inputtddft, err = 200, iostat = ios )
  ENDIF

  ! check input
  if (max_seconds < 0.1d0) call errore ('tddft_readin', ' wrong max_seconds', 1)
200 call errore('tddft_readin', 'reading inputtddft namelist', abs(ios))

  select case (verbosity)
     case('low')
       iverbosity = 1
     case('medium')
       iverbosity = 11
     case('high')
       iverbosity = 21
     case default
       call errore('tdddft_readin', 'verbosity can be ''low'', ''medium'' or ''high''', 1)
  end select

  ! convert to atomic units
  e_strength = e_strength * bohr_radius_angs ! change from 1/Ang to 1/Bohr
  dt = dt * 1.d-18 / (2.d0*au_sec)           ! change from femto-second to a.u. (Rydberg unit)
400 continue

#ifdef __MPI
  ! broadcast input variables  
  call tddft_bcast_input
#endif

END SUBROUTINE qepy_tddft_readin


!#ifdef __MPI
!!-----------------------------------------------------------------------
!SUBROUTINE tddft_bcast_input
!  !-----------------------------------------------------------------------
!  !
!  ! ... Broadcast input data to all processors 
!  !
!  USE mp_world,      ONLY : world_comm
!  USE mp,            ONLY : mp_bcast
!  USE io_files,      ONLY : prefix, tmp_dir
!  USE tddft_module

!  implicit none
!  integer, parameter :: root = 0    

!  call mp_bcast(job, root, world_comm)
!  call mp_bcast(prefix, root, world_comm)
!  call mp_bcast(tmp_dir, root, world_comm)
!  call mp_bcast(dt, root, world_comm)
!  call mp_bcast(e_strength, root, world_comm)
!  call mp_bcast(e_direction, root, world_comm)
!  call mp_bcast(conv_threshold, root, world_comm)
!  call mp_bcast(nstep, root, world_comm)
!  call mp_bcast(nupdate_Dnm , root, world_comm)
!  call mp_bcast(l_circular_dichroism, root, world_comm)
!  call mp_bcast(l_tddft_restart, root, world_comm)
!  call mp_bcast(iverbosity, root, world_comm)
!  call mp_bcast(max_seconds, root, world_comm)
!  call mp_bcast(molecule, root, world_comm)
!  call mp_bcast(ehrenfest, root, world_comm)
!  call mp_bcast(isave_rho, root, world_comm)
!  call mp_bcast(wavepacket, root, world_comm)
!  call mp_bcast(wp_pos, root, world_comm)
!  call mp_bcast(wp_d, root, world_comm)
!  call mp_bcast(wp_ekin, root, world_comm)

!END SUBROUTINE tddft_bcast_input
!#endif
  

!!-----------------------------------------------------------------------
SUBROUTINE qepy_tddft_allocate
  !-----------------------------------------------------------------------
  !
  ! ... Allocate memory for TDDFT
  !
  USE tddft_module
  USE klist,         ONLY : nkstot
  USE wvfct,         ONLY : btype, nbndx
  USE tddft_module

  implicit none

  ! needed by sum_band
  IF (.NOT.ALLOCATED(btype)) THEN
  allocate(btype(nbndx,nkstot))
  btype = 1
  ENDIF
    
END SUBROUTINE qepy_tddft_allocate


!!-----------------------------------------------------------------------
!SUBROUTINE tddft_summary
!  !-----------------------------------------------------------------------
!  !
!  ! ... Print a short summary of the calculation
!  !
!  USE io_global,        ONLY : stdout
!  USE lsda_mod,         ONLY : nspin
!  USE wvfct,            ONLY : nbnd
!  USE paw_variables,    ONLY : okpaw
!  USE uspp,             ONLY : okvan
!  USE tddft_module
!  implicit none
!  integer :: is
 
!  write(stdout,*)

!  write(stdout,'(5X,''Calculation type      : '',A12)') job
!  if (molecule) then
!     write(stdout,'(5X,''System is             : molecule'')')
!  else
!     write(stdout,'(5X,''System is             : crystal'')')
!  endif
!  if (ehrenfest) then
!     write(stdout,'(5X,''Ehrenfest dynamics'')')
!     if (okpaw .or. okvan) call infomsg('tddft_summary', 'Ehrenfest dynamics not yet supported with USPP and PAW')
!  endif
!  write(stdout,'(5X,''Number or steps       : '',I12)') nstep
!  write(stdout,'(5X,''Time step             : '',F12.4,'' rydberg_atomic_time'')') dt
!  write(stdout,'(5X,''Electric field dir.   : '',I12,'' (1=x,2=y,3=z)'')') e_direction
!  write(stdout,'(5X,''Electric field impulse: '',F12.4,'' bohrradius^-1'')') e_strength

!  write(stdout,*)

!  if (tfixed_occ) then
!     write(stdout,'(5X,''Occupations from input:'')')
!     do is = 1, nspin
!       write(stdout,'(5X,''ispin='',I1,'': '')',advance='no') is
!       write(stdout,'(10(F5.2,2X))') f_inp(1:nbnd,is)
!     enddo
!    write(stdout,*)
!  endif
     
!  flush( stdout )

!END SUBROUTINE tddft_summary
  
  

!!-----------------------------------------------------------------------
!SUBROUTINE tddft_openfil
!  !-----------------------------------------------------------------------
!  !
!  ! ... Open files needed for TDDFT
!  !
!  USE tddft_module   
!  USE wvfct,            ONLY : nbnd, npwx
!  USE ldaU,             ONLY : lda_plus_U, nwfcU
!  USE io_files,         ONLY : iunhub, iunwfc,nwordwfcU, nwordwfc, seqopn
!  USE noncollin_module, ONLY : npol
!  USE buffers,          ONLY : open_buffer
!  USE control_flags,    ONLY : io_level    
!  IMPLICIT NONE  
!  character*1, parameter :: dir(3) = (/'x', 'y', 'z'/)
!  logical :: exst

!  !
!  ! ... nwordwfc is the record length (IN REAL WORDS)
!  ! ... for the direct-access file containing wavefunctions
!  ! ... io_level > 0 : open a file; io_level <= 0 : open a buffer
!  !
!  nwordwfc = nbnd*npwx*npol
!  CALL open_buffer( iunwfc, 'wfc', nwordwfc, io_level, exst )

!  ! do not overwrite wfc
!  nwordwfc = nbnd*npwx*npol
!  CALL open_buffer( iunevcn, 'wfc'//dir(e_direction), nwordwfc, io_level, exst )

!  ! for restart
!  nwordtdwfc = nbnd*npwx*npol*2
!  CALL open_buffer( iuntdwfc, 'tmp'//dir(e_direction), nwordtdwfc, io_level, exst )

!  ! ... Needed for LDA+U
!  ! ... iunhub contains the (orthogonalized) atomic wfcs * S
!  nwordwfcU = npwx*nwfcU*npol
!  IF ( lda_plus_u ) &
!     CALL open_buffer( iunhub, 'hub', nwordwfcU, io_level, exst )

!END SUBROUTINE tddft_openfil


!-----------------------------------------------------------------------
SUBROUTINE qepy_tddft_closefil
  !-----------------------------------------------------------------------
  !
  ! ... Close files opened by TDDFT
  !
  USE ldaU,             ONLY : lda_plus_U  
  USE io_files,         ONLY : iunhub, iunwfc
  USE buffers,          ONLY : close_buffer
  USE tddft_module
  IMPLICIT NONE
  logical :: opnd

  call close_buffer( iunwfc, 'keep' )
  call close_buffer( iunevcn, 'keep' )
  call close_buffer( iuntdwfc, 'keep' )
  if ( lda_plus_u ) call close_buffer ( iunhub, status = 'keep' )
  inquire (unit = iunwfc, opened = opnd)
  inquire (unit = iunevcn, opened = opnd)
  inquire (unit = iuntdwfc, opened = opnd)

END SUBROUTINE qepy_tddft_closefil



!!-----------------------------------------------------------------------
!SUBROUTINE print_clock_tddft
!  !-----------------------------------------------------------------------
!  !
!  ! ... Print clocks
!  !
!  USE io_global,  ONLY : stdout
!  IMPLICIT NONE

!  write(stdout,*) '    Initialization:'
!  call print_clock ('tddft_setup')
!  write(stdout,*)
!  write(stdout,*) '    Linear response'
!  call print_clock ('greenf')
!  call print_clock ('cgsolve')
!  call print_clock ('ch_psi')
!  call print_clock ('h_psi')
!  call print_clock ('s_psi')
!  write(stdout,*) '    Real time evolution'
!  call print_clock ('updateH')
!  call print_clock ('dipole')
!  call print_clock ('quadrupole')
!  call print_clock ('circular')
!  write(stdout,*)
!  write(stdout,*) '    General routines'
!  call print_clock ('calbec')
!  call print_clock ('fft')
!  call print_clock ('ffts')
!  call print_clock ('fftw')
!  call print_clock ('cinterpolate')
!  call print_clock ('davcio')
!  call print_clock ('write_rec')
!  write(stdout,*)

!#ifdef __MPI
!  write(stdout,*) '    Parallel routines'
!  call print_clock ('reduce')  
!  call print_clock( 'fft_scatter' )
!  call print_clock( 'ALLTOALL' )
!  write(stdout,*)
!#endif
!  call print_clock ('TDDFT') 

!END SUBROUTINE print_clock_tddft



!!-----------------------------------------------------------------------
!SUBROUTINE tddft_memory_report
!  !-----------------------------------------------------------------------
!  !
!  ! ... Print estimated memory usage
!  !
!  USE io_global,                 ONLY : stdout
!  USE noncollin_module,          ONLY : npol
!  USE uspp,                      ONLY : nkb
!  USE fft_base,                  ONLY : dffts
!  USE pwcom
!  IMPLICIT NONE
!  integer, parameter :: Mb=1024*1024, complex_size=16, real_size=8

!  ! the conversions to double prevent integer overflow in very large run
!  write(stdout,'(5x,"Largest allocated arrays",5x,"est. size (Mb)",5x,"dimensions")')

!  write(stdout,'(8x,"KS wavefunctions at k     ",f10.2," Mb",5x,"(",i8,",",i5,")")') &
!     complex_size*nbnd*npol*DBLE(npwx)/Mb, npwx*npol,nbnd

!  write(stdout,'(8x,"First-order wavefunctions ",f10.2," Mb",5x,"(",i8,",",i5,",",i3")")') &
!     complex_size*nbnd*npol*DBLE(npwx)*10/Mb, npwx*npol,nbnd,10

!  write(stdout,'(8x,"Charge/spin density       ",f10.2," Mb",5x,"(",i8,",",i5,")")') &
!     real_size*dble(dffts%nnr)*nspin/Mb, dffts%nnr, nspin
  
!  write(stdout,'(8x,"NL pseudopotentials       ",f10.2," Mb",5x,"(",i8,",",i5,")")') &
!     complex_size*nkb*DBLE(npwx)/Mb, npwx, nkb
!  write(stdout,*)

!END SUBROUTINE tddft_memory_report


