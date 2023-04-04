!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE qepy_stop_run( exit_status, print_flag, what, finalize )
  !----------------------------------------------------------------------------
  !! Close all files and synchronize processes before stopping:
  !
  !! * exit_status = 0: successfull execution, remove temporary files;
  !! * exit_status =-1: code stopped by user request;
  !! * exit_status = 1: convergence not achieved.
  !
  !! Do not remove temporary files needed for restart.
  !
  !qepy --> 
  ! Also add some from pwscf and run_pwscf
  ! Merge and modify the mp_global.mp_global_end
  !qepy <-- 
  !
  USE io_global,          ONLY : stdout, ionode
  USE mp_global,          ONLY : mp_global_end
  USE environment,        ONLY : environment_end
  USE io_files,           ONLY : iuntmp, seqopn
  USE qmmm,               ONLY : qmmm_shutdown
  USE qexsd_module,       ONLY : qexsd_set_status
  USE mp,                 ONLY : mp_comm_free, mp_barrier, mp_start, mp_end, mp_stop, mp_count_nodes
  USE mp_world,           ONLY : world_comm
  USE mp_bands,           ONLY : inter_bgrp_comm, intra_bgrp_comm
  USE mp_pools,           ONLY : inter_pool_comm, intra_pool_comm
  USE mp_exx,             ONLY : inter_egrp_comm, intra_egrp_comm
  USE mp_images,          ONLY : inter_image_comm, intra_image_comm
  USE mp_orthopools,      ONLY : mp_stop_orthopools
  USE mp_diag,            ONLY : ortho_comm
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: exit_status
  INTEGER, INTENT(IN), OPTIONAL :: print_flag
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: what
  LOGICAL, OPTIONAL   :: finalize
  CHARACTER(LEN=256)  :: what_ ='config-nowf'
  INTEGER             :: iprint = 0
  LOGICAL             :: exst, opnd, lflag
#if defined(__MPI)
  INTEGER :: ierr
#endif
  !
  !qepy --> pwscf and run_pwscf
  CALL qexsd_set_status( exit_status )
  CALL qmmm_shutdown()
  if ( ortho_comm  /= 0 .and. ortho_comm  /= world_comm ) CALL laxlib_free_ortho_group()
  !qepy <-- pwscf and run_pwscf
  !
  IF ( PRESENT(what)) THEN
     IF (len_trim(what)>1) what_= trim(what)
  ENDIF

  IF (TRIM(what_) == 'no') THEN 
     WRITE( UNIT = stdout, FMT = '(/,5X,"Not output data")' )
  ELSE
     CALL punch( what_)
  ENDIF

  IF ( PRESENT(print_flag)) THEN
     iprint = print_flag
  ENDIF


  lflag = ( exit_status == 0 ) 
  IF ( lflag ) THEN
     ! 
     ! ... remove files needed only to restart
     !
     CALL seqopn( iuntmp, 'restart', 'UNFORMATTED', exst )
     CLOSE( UNIT = iuntmp, STATUS = 'DELETE' )
     !
     IF ( ionode ) THEN
        CALL seqopn( iuntmp, 'update', 'FORMATTED', exst )
        CLOSE( UNIT = iuntmp, STATUS = 'DELETE' )
        CALL seqopn( iuntmp, 'para', 'FORMATTED', exst )
        CLOSE( UNIT = iuntmp, STATUS = 'DELETE' )
     ENDIF
     !
  ENDIF
  !
  CALL close_files( lflag )
  !
  IF ( iprint > 0 .and. iprint<10 ) THEN
     CALL print_clock_pw()
  ENDIF
  !
  CALL clean_pw( .TRUE. )
  !
  IF ( iprint > 0 .and. iprint<10 ) THEN
     CALL environment_end( 'PWSCF' )
  ENDIF
  !
  !CALL mp_global_end()
  !-----------------------------------------------------------------------
  !qepy --> add mp_global_end
  if ( intra_egrp_comm  /= 0 .and. intra_egrp_comm  /= world_comm ) CALL mp_comm_free ( intra_egrp_comm )
  if ( inter_egrp_comm  /= 0 .and. inter_egrp_comm  /= world_comm ) CALL mp_comm_free ( inter_egrp_comm )
  if ( intra_bgrp_comm  /= 0 .and. intra_bgrp_comm  /= world_comm ) CALL mp_comm_free ( intra_bgrp_comm )
  if ( inter_bgrp_comm  /= 0 .and. inter_bgrp_comm  /= world_comm ) CALL mp_comm_free ( inter_bgrp_comm )
  if ( intra_pool_comm  /= 0 .and. intra_pool_comm  /= world_comm ) CALL mp_comm_free ( intra_pool_comm )
  if ( inter_pool_comm  /= 0 .and. inter_pool_comm  /= world_comm ) CALL mp_comm_free ( inter_pool_comm )
  if ( intra_image_comm /= 0 .and. intra_image_comm /= world_comm ) CALL mp_comm_free ( intra_image_comm )
  if ( inter_image_comm /= 0 .and. inter_image_comm /= world_comm ) CALL mp_comm_free ( inter_image_comm )
  CALL mp_stop_orthopools( ) ! cleans orthopools if used in exx
  CALL mp_barrier( world_comm )
  CALL mp_end ( world_comm )
#if defined(__MPI)
  IF (present(finalize)) THEN
     IF (finalize) THEN
        CALL mpi_finalize(ierr)
        IF (ierr/=0) CALL mp_stop( 8002 )
     ENDIF
  END IF
#endif
  !qepy <-- add mp_global_end
  !-----------------------------------------------------------------------
  !
END SUBROUTINE qepy_stop_run
!
!-----------------------------------------
!SUBROUTINE do_stop( exit_status )
  !!---------------------------------------
  !!! Stop the run.
  !!
  !IMPLICIT NONE
  !!
  !INTEGER, INTENT(IN) :: exit_status
  !!
  !IF ( exit_status == -1 ) THEN
     !! -1 is not an acceptable value for stop in fortran;
     !! convert it to 255
     !STOP 255
  !ELSEIF ( exit_status == 0 ) THEN
     !STOP
  !ELSEIF ( exit_status == 1 ) THEN
     !STOP 1
  !ELSEIF ( exit_status == 2 ) THEN
     !STOP 2
  !ELSEIF ( exit_status == 3 ) THEN
     !STOP 3
  !ELSEIF ( exit_status == 4 ) THEN
     !STOP 4
  !ELSEIF ( exit_status == 255 ) THEN
     !STOP 255
  !ELSEIF ( exit_status == 254 ) THEN
     !STOP 254
  !ELSE
     !! unimplemented value
     !STOP 128
  !ENDIF
  !!
!END SUBROUTINE do_stop
!!
!!----------------------------------------------------------------------------
!SUBROUTINE closefile()
  !!----------------------------------------------------------------------------
  !!! Close all files and synchronize processes before stopping.  
  !!! Called by "sigcatch" when it receives a signal.
  !!
  !USE io_global,  ONLY :  stdout
  !!
  !WRITE( stdout,'(5X,"Signal Received, stopping ... ")')
  !!
  !CALL stop_run( 255 )
  !!
  !RETURN
  !!
!END SUBROUTINE closefile
