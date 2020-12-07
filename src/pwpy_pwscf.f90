!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE pwpy_pwscf(infile, my_world_comm)
  !! Author: Paolo Giannozzi
  !
  !! Version: v6.1
  !
  !! License: GNU
  !
  !! Summary: Main program calling one or more instances of Plane Wave Self-Consistent Field code
  !
  !! This is the main program for executable "pw.x".  |
  !! * If called as "pw.x -ipi server-address" or "pw.x --ipi server-address",
  !! works in "server" mode, calls [[run_driver]].
  !! * If called as "manypw.x" via a link, works in "manypw" mode, runs many
  !! instances (images) of pw.x (see [[run_manypw]])
  !! * If called as "dist.x" via a link, works in "dry run" mode, computes
  !! distances, angles, neighbors, writes to file "dist.out" and stops. 
  !! Otherwise: see [[run_pwscf]]
  !!
  !! @Note
  !! 10/01/17 Samuel Ponce: Add Ford documentation
  !! @endnote
  !!
  !! @warning
  !! Example of Warning
  !!
  !! @todo
  !! Have automatic parallelisation. 
  !!
  !! @bug
  !! No bug.
  !!
  USE environment,          ONLY : environment_start
  USE mp_global,            ONLY : mp_startup
  USE mp_world,             ONLY : world_comm
  USE mp_pools,             ONLY : intra_pool_comm
  USE mp_bands,             ONLY : intra_bgrp_comm, inter_bgrp_comm
  USE mp_diag,              ONLY : mp_start_diag
  USE mp_exx,               ONLY : negrp
  USE read_input,           ONLY : read_input_file
  USE command_line_options, ONLY : input_file_, command_line, ndiag_
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*) :: infile
  INTEGER, INTENT(IN), OPTIONAL :: my_world_comm
  !
  CHARACTER(len=256) :: srvaddress
  !! Get the address of the server 
  CHARACTER(len=256) :: get_server_address
  !! Get the address of the server 
  INTEGER :: exit_status
  !! Status at exit
  LOGICAL :: use_images, do_diag_in_band_group = .TRUE.
  !! true if running "manypw.x"
  LOGICAL, EXTERNAL :: matches
  !! checks if first string is contained in the second
  !
  IF ( PRESENT(my_world_comm)) THEN
     CALL mp_startup(my_world_comm=my_world_comm, start_images=.TRUE. )
  ELSE
     CALL mp_startup( start_images=.TRUE. )
  ENDIF
  !
  IF( negrp > 1 .OR. do_diag_in_band_group ) THEN
     ! used to be the default : one diag group per bgrp
     ! with strict hierarchy: POOL > BAND > DIAG
     ! if using exx groups from mp_exx still use this diag method
     CALL mp_start_diag( ndiag_, world_comm, intra_bgrp_comm, &
                         do_distr_diag_inside_bgrp_ = .TRUE. )
  ELSE
     ! new default: one diag group per pool ( individual k-point level )
     ! with band group and diag group both being children of POOL comm
     CALL mp_start_diag( ndiag_, world_comm, intra_pool_comm, &
                         do_distr_diag_inside_bgrp_ = .FALSE. )
  ENDIF
  !
  CALL set_mpi_comm_4_solvers( intra_pool_comm, intra_bgrp_comm, &
                               inter_bgrp_comm )
  !
  CALL environment_start( 'PWSCF' )
  !
  ! ... Check if running standalone or in "driver" mode
  !
  !srvaddress = get_server_address( command_line ) 
  !
  ! ... Check if running standalone or in "manypw" mode
  !
  !use_images = matches( 'manypw.x', command_line )
  !
  ! ... Perform actual calculation
  !
  !IF ( TRIM(srvaddress) == ' ' ) THEN
    !! When running standalone:
    !IF ( use_images ) THEN
       !! as manypw.x
       !CALL run_manypw( )
       !CALL run_pwscf( exit_status )
       !!
     !ELSE
       !! as pw.x
       !CALL read_input_file( 'PW', input_file_ )
       !CALL run_pwscf( exit_status )
       !!
    !ENDIF
  !ELSE
     !! When running as library
     !!
     !CALL read_input_file('PW+iPi', input_file_ )
     !CALL run_driver( srvaddress, exit_status )
     !!
  !ENDIF
  !
  input_file_=trim(infile)
  CALL read_input_file( 'PW', input_file_ )
  call pwpy_run_pwscf(exit_status)
END SUBROUTINE pwpy_pwscf
   !
SUBROUTINE pwpy_pwscf_finalise()
   IMPLICIT NONE
   INTEGER :: exit_status

   CALL laxlib_free_ortho_group()
   CALL pwpy_stop_run( exit_status )
   !CALL do_stop( exit_status )
END SUBROUTINE pwpy_pwscf_finalise
