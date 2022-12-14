!
! Copyright (C) 2005-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE oldxml_xml_io_base
  !----------------------------------------------------------------------------
  !
  ! ... this module contains some common subroutines used to read and write
  ! ... in XML format the data produced by Quantum ESPRESSO package
  !
  ! ... written by Carlo Sbraccia (2005)
  ! ... modified by Andrea Ferretti (2006-08)
  !
  USE iotk_module
  !
  USE kinds,     ONLY : DP
  !USE io_files,  ONLY : tmp_dir, prefix, iunpun, xmlpun, check_file_exist, &
  USE io_files,  ONLY : tmp_dir, prefix, iunpun, check_file_exist, &
       create_directory
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp,        ONLY : mp_bcast, mp_sum, mp_get, mp_put, mp_max, mp_rank, &
       mp_size
  USE parser,    ONLY : version_compare
  !
  IMPLICIT NONE
  PRIVATE
  !
  CHARACTER(iotk_attlenx)  :: attr
  LOGICAL,       SAVE      :: rho_binary = .TRUE.
  !
  CHARACTER (LEN=13), PARAMETER :: xmlpun      = 'data-file.xml'
  CHARACTER(len=256) :: qexml_version = ' '       ! the format of the current qexml datafile 
  LOGICAL            :: qexml_version_init = .FALSE.  ! whether the fmt has been read or not
  !
  PUBLIC :: xmlpun, qexml_version, qexml_version_init
  PUBLIC :: rho_binary
  PUBLIC :: attr
  !
  PUBLIC :: read_wfc, write_wfc, read_rho, write_rho, &
       save_print_counter, read_print_counter, restart_dir
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    FUNCTION restart_dir( outdir, runit )
      !------------------------------------------------------------------------
      !
      ! KNK_nimage
      ! USE mp_images, ONLY:  my_image_id
      CHARACTER(LEN=256)           :: restart_dir
      CHARACTER(LEN=*), INTENT(IN) :: outdir
      INTEGER,          INTENT(IN) :: runit
      !
      CHARACTER(LEN=256)         :: dirname
      INTEGER                    :: strlen
      CHARACTER(LEN=6), EXTERNAL :: int_to_char
      !
      ! ... main restart directory
      !
      dirname = TRIM( prefix ) // '_' // TRIM( int_to_char( runit ) )// '.save/'
      !
      IF ( LEN( outdir ) > 1 ) THEN
         !
         strlen = INDEX( outdir, ' ' ) - 1
         !
         dirname = outdir(1:strlen) // '/' // dirname
         !
      END IF
      !
      restart_dir = TRIM( dirname )
      !
      RETURN
      !
    END FUNCTION restart_dir
    !
    !------------------------------------------------------------------------
    FUNCTION check_restartfile( outdir, ndr )
      !------------------------------------------------------------------------
      !
      USE mp_images, ONLY : intra_image_comm
      !
      IMPLICIT NONE
      !
      LOGICAL                      :: check_restartfile
      INTEGER,          INTENT(IN) :: ndr
      CHARACTER(LEN=*), INTENT(IN) :: outdir
      CHARACTER(LEN=256)           :: filename
      LOGICAL                      :: lval
      !
      !
      filename = restart_dir( outdir, ndr )
      !
      IF ( ionode ) THEN
         !
         filename = TRIM( filename ) // '/' // TRIM( xmlpun )
         !
         INQUIRE( FILE = TRIM( filename ), EXIST = lval )
         !
      END IF
      !
      CALL mp_bcast( lval, ionode_id, intra_image_comm )
      !
      check_restartfile = lval
      !
      RETURN
      !
    END FUNCTION check_restartfile
    !
    !------------------------------------------------------------------------
    SUBROUTINE save_print_counter( iter, outdir, wunit )
      !------------------------------------------------------------------------
      !
      ! ... a counter indicating the last successful printout iteration is saved
      !
      USE mp_images, ONLY : intra_image_comm
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(IN) :: iter
      CHARACTER(LEN=*), INTENT(IN) :: outdir
      INTEGER,          INTENT(IN) :: wunit
      !
      INTEGER            :: ierr
      CHARACTER(LEN=256) :: filename, dirname
      !
      !
      dirname = restart_dir( outdir, wunit )
      !
      CALL create_directory( TRIM( dirname ) )
      !
      IF ( ionode ) THEN
         !
         filename = TRIM( dirname ) // 'print_counter.xml'
         !
         CALL iotk_open_write( iunpun, FILE = filename, &
                             & ROOT = "PRINT_COUNTER",  IERR = ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'save_print_counter', &
                   'cannot open restart file for writing', ierr )
      !
      IF ( ionode ) THEN
         !
         CALL iotk_write_begin( iunpun, "LAST_SUCCESSFUL_PRINTOUT" )
         CALL iotk_write_dat(   iunpun, "STEP", iter )
         CALL iotk_write_end(   iunpun, "LAST_SUCCESSFUL_PRINTOUT" )
         !
         CALL iotk_close_write( iunpun )
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE save_print_counter
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_print_counter( nprint_nfi, outdir, runit )
      !------------------------------------------------------------------------
      !
      ! ... the counter indicating the last successful printout iteration 
      ! ... is read here
      !
      USE mp_images, ONLY : intra_image_comm
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(OUT) :: nprint_nfi
      CHARACTER(LEN=*), INTENT(IN)  :: outdir
      INTEGER,          INTENT(IN)  :: runit
      !
      INTEGER            :: ierr
      CHARACTER(LEN=256) :: filename, dirname
      !
      !
      dirname = restart_dir( outdir, runit )
      !
      IF ( ionode ) THEN
         !
         filename = TRIM( dirname ) // 'print_counter.xml'
         !
         CALL iotk_open_read( iunpun, FILE = filename, IERR = ierr )
         !
         IF ( ierr > 0 ) THEN
            !
            nprint_nfi = -1
            !
         ELSE
            !
            CALL iotk_scan_begin( iunpun, "LAST_SUCCESSFUL_PRINTOUT" )
            CALL iotk_scan_dat(   iunpun, "STEP", nprint_nfi )
            CALL iotk_scan_end(   iunpun, "LAST_SUCCESSFUL_PRINTOUT" )
            !
            CALL iotk_close_read( iunpun )
            !
         END IF
         !
      END IF
      !
      CALL mp_bcast( nprint_nfi, ionode_id, intra_image_comm )
      !
      RETURN
      !
    END SUBROUTINE read_print_counter   
    !
    !------------------------------------------------------------------------
    SUBROUTINE set_kpoints_vars( ik, nk, kunit, ngwl, igl, &
                                 ngroup, ikt, iks, ike, igwx, ipmask, ipsour, &
                                 ionode, root_in_group, intra_group_comm, inter_group_comm, parent_group_comm )
      !------------------------------------------------------------------------
      !
      ! ... set working variables for k-point index (ikt) and 
      ! ... k-points number (nkt)
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN)  :: ik, nk, kunit
      INTEGER, INTENT(IN)  :: ngwl, igl(:)
      INTEGER, INTENT(OUT) :: ngroup
      INTEGER, INTENT(OUT) :: ikt, iks, ike, igwx
      INTEGER, INTENT(OUT) :: ipmask(:), ipsour
      LOGICAL, INTENT(IN)  :: ionode
      INTEGER, INTENT(IN)  :: root_in_group, intra_group_comm, inter_group_comm, parent_group_comm
      !
      INTEGER :: ierr, i
      INTEGER :: nkl, nkr, nkbl, nkt
      INTEGER :: nproc_parent, nproc_group, my_group_id, me_in_group, me_in_parent, io_in_parent
      !
      nproc_parent = mp_size( parent_group_comm )
      nproc_group  = mp_size( intra_group_comm )
      my_group_id  = mp_rank( inter_group_comm )
      me_in_group  = mp_rank( intra_group_comm )
      me_in_parent = mp_rank( parent_group_comm )
      !
      ! find the ID (io_in_parent) of the io PE ( where ionode == .true. )
      !
      io_in_parent = 0
      IF( ionode ) io_in_parent = me_in_parent
      CALL mp_sum( io_in_parent, parent_group_comm )
      !
      ikt = ik
      nkt = nk
      !
      ! ... find out the number of pools
      !
      ngroup = nproc_parent / nproc_group 
      !
      ! ... find out number of k points blocks
      !
      nkbl = nkt / kunit  
      !
      ! ... k points per pool
      !
      nkl = kunit * ( nkbl / ngroup )
      !
      ! ... find out the reminder
      !
      nkr = ( nkt - nkl * ngroup ) / kunit
      !
      ! ... Assign the reminder to the first nkr pools
      !
      IF ( my_group_id < nkr ) nkl = nkl + kunit
      !
      ! ... find out the index of the first k point in this pool
      !
      iks = nkl * my_group_id + 1
      !
      IF ( my_group_id >= nkr ) iks = iks + nkr * kunit
      !
      ! ... find out the index of the last k point in this pool
      !
      ike = iks + nkl - 1
      !
      ipmask = 0
      ipsour = io_in_parent
      !
      ! ... find out the index of the processor which collect the data 
      ! ... in the pool of ik
      !
      IF ( ngroup > 1 ) THEN
         !
         IF ( ( ikt >= iks ) .AND. ( ikt <= ike ) ) THEN
            !
            IF ( me_in_group == root_in_group ) ipmask( me_in_parent + 1 ) = 1
            !
         END IF
         !
         ! ... Collect the mask for all proc in the image
         !
         CALL mp_sum( ipmask, parent_group_comm )
         !
         DO i = 1, nproc_parent
            !
            IF( ipmask(i) == 1 ) ipsour = ( i - 1 )
            !
         END DO
         !
      END IF
      !
      igwx = 0
      ierr = 0
      !
      IF ( ( ikt >= iks ) .AND. ( ikt <= ike ) ) THEN
         !
         IF ( ngwl > SIZE( igl ) ) THEN
            !
            ierr = 1
            !
         ELSE
            !
            igwx = MAXVAL( igl(1:ngwl) )
            !
         END IF
         !
      END IF
      !
      ! ... get the maximum index within the pool
      !
      CALL mp_max( igwx, intra_group_comm )
      !
      ! ... now notify all procs if an error has been found 
      !
      CALL mp_max( ierr, parent_group_comm )
      !
      CALL errore( 'set_kpoint_vars ', 'wrong size ngl', ierr )
      !
      IF ( ipsour /= io_in_parent ) &
         CALL mp_get( igwx, igwx, me_in_parent, io_in_parent, ipsour, 1, parent_group_comm )
      !
      RETURN
      !
    END SUBROUTINE set_kpoints_vars
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_rho( dirname, rho, nspin, extension )
      !------------------------------------------------------------------------
      !
      ! ... this routine writes the charge-density in xml format into the
      ! ... $dirname directory - $dirname must exist and end with '/'
      !
      USE fft_base, ONLY : dfftp
      USE mp_bands, ONLY : intra_bgrp_comm, inter_bgrp_comm
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(IN)           :: nspin
      REAL(DP),         INTENT(IN)           :: rho(dfftp%nnr,nspin)
      CHARACTER(LEN=*), INTENT(IN)           :: dirname
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extension
      !
      CHARACTER(LEN=256)    :: file_base
      CHARACTER(LEN=6)      :: ext
      REAL(DP), ALLOCATABLE :: rhoaux(:)
      !
      !
      ext = ' '
      IF ( PRESENT( extension ) ) ext = '.' // TRIM( extension )
      !
      file_base = TRIM( dirname ) // 'charge-density' // TRIM( ext )
      !
      IF ( nspin == 1 ) THEN
         !
         CALL write_rho_xml( file_base, rho(:,1), dfftp, ionode, inter_bgrp_comm )
         !
      ELSE IF ( nspin == 2 ) THEN
         !
         ALLOCATE( rhoaux( dfftp%nnr ) )
         !
         rhoaux(:) = rho(:,1) + rho(:,2)
         !
         CALL write_rho_xml( file_base, rhoaux, dfftp, ionode, inter_bgrp_comm )
         !
         file_base = TRIM( dirname ) // 'spin-polarization' // TRIM( ext )
         !
         rhoaux(:) = rho(:,1) - rho(:,2)
         !
         CALL write_rho_xml( file_base, rhoaux,  dfftp, ionode, inter_bgrp_comm )
         !
         DEALLOCATE( rhoaux )
         !
      ELSE IF ( nspin == 4 ) THEN
         !
         CALL write_rho_xml( file_base, rho(:,1), dfftp, ionode, inter_bgrp_comm )
         !
         file_base = TRIM( dirname ) // 'magnetization.x' // TRIM( ext )
         !
         CALL write_rho_xml( file_base, rho(:,2), dfftp, ionode, inter_bgrp_comm )
         !
         file_base = TRIM( dirname ) // 'magnetization.y' // TRIM( ext )
         !
         CALL write_rho_xml( file_base, rho(:,3), dfftp, ionode, inter_bgrp_comm )
         !
         file_base = TRIM( dirname ) // 'magnetization.z' // TRIM( ext )
         !
         CALL write_rho_xml( file_base, rho(:,4), dfftp, ionode, inter_bgrp_comm )
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE write_rho
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_rho( dirname, rho, nspin, extension)
      !------------------------------------------------------------------------
      !
      ! ... this routine reads the charge-density in xml format from the
      ! ... files saved into the '.save' directory
      !
      USE fft_base,  ONLY : dfftp
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(IN)           :: nspin
      CHARACTER(LEN=*), INTENT(IN)           :: dirname
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extension
      REAL(DP),         INTENT(OUT)          :: rho(dfftp%nnr,nspin)
      !
      CHARACTER(LEN=256)  :: file_base
      CHARACTER(LEN=6)    :: ext
      REAL(DP), ALLOCATABLE :: rhoaux(:)
      !
      ext = ' '
      IF ( PRESENT( extension ) ) ext = '.' // TRIM( extension )
      !
      file_base = TRIM( dirname ) // 'charge-density' // TRIM( ext )
      CALL read_rho_xml ( file_base, dfftp, rho(:,1) ) 
      !
      IF ( nspin == 2 ) THEN
         !
         rho(:,2) = rho(:,1)
         !
         ALLOCATE( rhoaux( dfftp%nnr ) )
         !
         file_base = TRIM( dirname ) // 'spin-polarization' // TRIM( ext )
         CALL read_rho_xml ( file_base, dfftp, rhoaux ) 
         !
         rho(:,1) = 0.5D0*( rho(:,1) + rhoaux(:) )
         rho(:,2) = 0.5D0*( rho(:,2) - rhoaux(:) )
         !
         DEALLOCATE( rhoaux )
         !
      ELSE IF ( nspin == 4 ) THEN
         !
         file_base = TRIM( dirname ) // 'magnetization.x' // TRIM( ext )
         CALL read_rho_xml ( file_base, dfftp, rho(:,2) ) 
         !
         file_base = TRIM( dirname ) // 'magnetization.y' // TRIM( ext )
         CALL read_rho_xml ( file_base, dfftp, rho(:,3) ) 
         !
         file_base = TRIM( dirname ) // 'magnetization.z' // TRIM( ext )
         CALL read_rho_xml ( file_base, dfftp, rho(:,4) ) 
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE read_rho
    !------------------------------------------------------------------------
    SUBROUTINE write_rho_xml( rho_file_base, rho, fft_desc, ionode, inter_group_comm )
      !------------------------------------------------------------------------
      !
      ! ... Writes charge density rho, one plane at a time.
      ! ... If ipp and npp are specified, planes are collected one by one from
      ! ... all processors, avoiding an overall collect of the charge density
      ! ... on a single proc.
      !
#if defined __HDF5
      USE hdf5_qe,  ONLY  : write_rho_hdf5, h5fclose_f, &
                            prepare_for_writing_final, add_attributes_hdf5, rho_hdf5_write  
#endif
      USE fft_types
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*),  INTENT(IN) :: rho_file_base
      REAL(DP),          INTENT(IN) :: rho(:)
      TYPE(fft_type_descriptor),INTENT(IN) :: fft_desc
      INTEGER,           INTENT(IN) :: inter_group_comm
      LOGICAL,           INTENT(IN) :: ionode
      !
      INTEGER               :: nr1,nr2,nr3, nr1x, nr2x,nr3x
      INTEGER               :: rhounit, ierr, i, j, jj, k, kk, ldr, ip
      CHARACTER(LEN=256)    :: rho_file
      CHARACTER(LEN=256)    :: rho_file_hdf5
      CHARACTER(LEN=10)     :: rho_extension
      REAL(DP), ALLOCATABLE :: rho_plane(:)
      INTEGER,  ALLOCATABLE :: kowner(:)
      INTEGER               :: my_group_id, me_group, me_group2, me_group3, &
                                            nproc_group, nproc_group2, nproc_group3, &
                                            io_group_id, io_group2, io_group3
      INTEGER,  EXTERNAL    :: find_free_unit
      !

      my_group_id = mp_rank( inter_group_comm )

      me_group = fft_desc%mype ; me_group2 = fft_desc%mype2 ; me_group3 = fft_desc%mype3
      nproc_group = fft_desc%nproc ; nproc_group2 = fft_desc%nproc2 ; nproc_group3 = fft_desc%nproc3
      !
      nr1  = fft_desc%nr1  ; nr2  = fft_desc%nr2  ; nr3  = fft_desc%nr3
      nr1x = fft_desc%nr1x ; nr2x = fft_desc%nr2x ; nr3x = fft_desc%nr3x
      !
      rho_extension = '.dat'
      IF ( .NOT. rho_binary ) rho_extension = '.xml'
      !
      rho_file = TRIM( rho_file_base ) // TRIM( rho_extension )
      rhounit = find_free_unit ()
      !
      IF ( ionode ) THEN 
#if defined  __HDF5
         rho_file_hdf5 = TRIM( rho_file_base ) // '.hdf5'
         CALL prepare_for_writing_final(rho_hdf5_write, 0 ,rho_file_hdf5)
         CALL add_attributes_hdf5(rho_hdf5_write,nr1,"nr1")
         CALL add_attributes_hdf5(rho_hdf5_write,nr2,"nr2")
         CALL add_attributes_hdf5(rho_hdf5_write,nr3,"nr3")
#else
         CALL iotk_open_write( rhounit, FILE = rho_file,  BINARY = rho_binary, IERR = ierr )
         CALL errore( 'write_rho_xml', 'cannot open ' // TRIM( rho_file ) // ' file for writing', ierr )
#endif
      END IF 
      !
#if !defined __HDF5
      IF ( ionode ) THEN
         !
         CALL iotk_write_begin( rhounit, "CHARGE-DENSITY" )
         !
         CALL iotk_write_attr( attr, "nr1", nr1, FIRST = .TRUE. )
         CALL iotk_write_attr( attr, "nr2", nr2 )
         CALL iotk_write_attr( attr, "nr3", nr3 )
         !
         CALL iotk_write_empty( rhounit, "INFO", attr )
         !
      END IF
#endif
      !
      ALLOCATE( rho_plane( nr1*nr2 ) )
      ALLOCATE( kowner( nr3 ) )
      !
      ! ... find the index of the group (pool) that will write rho
      !
      io_group_id = 0
      !
      IF ( ionode ) io_group_id = my_group_id
      !
      CALL mp_sum( io_group_id, fft_desc%comm )
      CALL mp_sum( io_group_id, inter_group_comm ) ! io_group_id is the (pool) group that contains the ionode
      !
      ! ... find the index of the ionode within Y and Z  groups
      !
      io_group2 = 0 ; IF ( ionode ) io_group2 = me_group2
      CALL mp_sum( io_group2, fft_desc%comm )  ! io_group2 is the group index of the ionode in the Y group (nproc2)
      io_group3 = 0 ; IF ( ionode ) io_group3 = me_group3
      CALL mp_sum( io_group3, fft_desc%comm )  ! io_group3 is the group index of the ionode in the Z group (nproc3)
      !
      ! ... find out the owner of each "z" plane
      !
      DO ip = 1, nproc_group3
         !
         kowner( (fft_desc%i0r3p(ip)+1):(fft_desc%i0r3p(ip)+fft_desc%nr3p(ip)) ) = ip - 1
         !
      END DO
      !
      ldr = nr1x*fft_desc%my_nr2p
      !
      IF ( ( my_group_id == io_group_id ) ) THEN ! only the group of ionode collects and writes the data
         !
         DO k = 1, nr3
            !
            !  Only one subgroup write the charge density
            ! 
            rho_plane = 0.d0
            IF( ( kowner(k) == me_group3 ) ) THEN
               !
               kk = k - fft_desc%my_i0r3p
               ! 
               DO jj = 1, fft_desc%my_nr2p
                  !
                  j = jj + fft_desc%my_i0r2p
                  DO i = 1, nr1
                     !
                     rho_plane(i+(j-1)*nr1) = rho(i+(jj-1)*nr1x+(kk-1)*ldr)
                     !
                  END DO
                  !
               END DO
               call mp_sum(rho_plane, fft_desc%comm2 ) ! collect the data over the Y group (nproc2)
               !
            END IF
            !
            ! if this processor is in the same comm3 group as ionode (me_group2==io_group2) 
            IF ( kowner(k) /= io_group3 .and. me_group2==io_group2) & 
               CALL mp_get( rho_plane, rho_plane, me_group3, io_group3, kowner(k), k, fft_desc%comm3 )
            !
            IF ( ionode ) THEN
#if defined __HDF5
            CALL write_rho_hdf5(rho_hdf5_write,k,rho_plane)
#else
               CALL iotk_write_dat( rhounit, "z" // iotk_index( k ), rho_plane )
#endif
            ENDIF
            !
         END DO
         !
      END IF
      !
      DEALLOCATE( rho_plane )
      DEALLOCATE( kowner )
      !
      IF ( ionode ) THEN
#if defined __HDF5
         CALL h5fclose_f(rho_hdf5_write%file_id,ierr)
#else
   !
         CALL iotk_write_end( rhounit, "CHARGE-DENSITY" )
         !
         CALL iotk_close_write( rhounit )
#endif       
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE write_rho_xml
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_rho_xml( rho_file_base, fft_desc, rho )
      !------------------------------------------------------------------------
      !
      ! ... Reads charge density rho, one plane at a time, to avoid 
      ! ... collecting the entire charge density on a single processor
      !
      USE mp_images, ONLY : intra_image_comm
#if defined __HDF5
      USE hdf5_qe,   ONLY : read_rho_hdf5, read_attributes_hdf5, &
           prepare_for_reading_final, h5fclose_f, rho_hdf5_write, hdf5_type
#endif
      USE fft_types
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*),  INTENT(IN)  :: rho_file_base
      TYPE(fft_type_descriptor),INTENT(IN) :: fft_desc
      REAL(DP),          INTENT(OUT) :: rho(:)
      !
      INTEGER               :: rhounit, ierr, i, j, jj, k, kk, ldr, ip
      INTEGER               :: nr( 3 ), nr1_, nr2_, nr3_, nr1, nr2, nr3, nr1x, nr2x, nr3x
      INTEGER               :: me_group, me_group2, me_group3, &
                               nproc_group, nproc_group2, nproc_group3
      CHARACTER(LEN=256)    :: rho_file
      CHARACTER(LEN=256)    :: rho_file_hdf5
      REAL(DP), ALLOCATABLE :: rho_plane(:)
      INTEGER,  ALLOCATABLE :: kowner(:)
      LOGICAL               :: exst
      INTEGER,  EXTERNAL    :: find_free_unit
#if defined(__HDF5)
      TYPE(hdf5_type),ALLOCATABLE   :: h5desc
#endif
      !
      me_group = fft_desc%mype ; me_group2 = fft_desc%mype2 ; me_group3 = fft_desc%mype3
      nproc_group = fft_desc%nproc ; nproc_group2 = fft_desc%nproc2 ; nproc_group3 = fft_desc%nproc3
      !
      nr1  = fft_desc%nr1  ; nr2  = fft_desc%nr2  ; nr3  = fft_desc%nr3
      nr1x = fft_desc%nr1x ; nr2x = fft_desc%nr2x ; nr3x = fft_desc%nr3x
      !
#if defined(__HDF5)
      rho_file_hdf5 = TRIM( rho_file_base ) // '.hdf5'
      exst = check_file_exist(TRIM(rho_file_hdf5))
      IF ( .NOT. exst ) CALL errore ('read_rho_xml', 'searching for '// TRIM(rho_file_hdf5),10)
#else 
      rhounit = find_free_unit ( )
      rho_file = TRIM( rho_file_base ) // ".dat"
      exst = check_file_exist( TRIM(rho_file) ) 
      !
      IF ( .NOT. exst ) CALL errore('read_rho_xml', 'searching for '//TRIM(rho_file), 10)
#endif
      !
      IF ( ionode ) THEN
#if defined (__HDF5)
         ALLOCATE ( h5desc)
         CALL prepare_for_reading_final(h5desc, 0 ,rho_file_hdf5)
         CALL read_attributes_hdf5(h5desc, nr1_,"nr1")
         CALL read_attributes_hdf5(h5desc, nr2_,"nr2")
         CALL read_attributes_hdf5(h5desc, nr3_,"nr3")
         nr = [nr1_,nr2_,nr3_]
#else
         CALL iotk_open_read( rhounit, FILE = rho_file, IERR = ierr )
         CALL errore( 'read_rho_xml', 'cannot open ' // TRIM( rho_file ) // ' file for reading', ierr )
         CALL iotk_scan_begin( rhounit, "CHARGE-DENSITY" )
         !
         CALL iotk_scan_empty( rhounit, "INFO", attr )
         !
         CALL iotk_scan_attr( attr, "nr1", nr(1) )
         CALL iotk_scan_attr( attr, "nr2", nr(2) )
         CALL iotk_scan_attr( attr, "nr3", nr(3) )
#endif
         !
         IF ( nr1 /= nr(1) .OR. nr2 /= nr(2) .OR. nr3 /= nr(3) ) &
            CALL errore( 'read_rho_xml', 'dimensions do not match', 1 )
         !
      END IF
      !
      ALLOCATE( rho_plane( nr1*nr2 ) )
      ALLOCATE( kowner( nr3 ) )
      !
      DO ip = 1, nproc_group3
         !
         kowner( (fft_desc%i0r3p(ip)+1):(fft_desc%i0r3p(ip)+fft_desc%nr3p(ip)) ) = ip - 1
         !
      END DO
      !
      ldr = nr1x*fft_desc%my_nr2p
      !
      ! ... explicit initialization to zero is needed because the physical
      ! ... dimensions of rho may exceed the true size of the FFT grid 
      !
      rho(:) = 0.0_DP
      !
      DO k = 1, nr3
         !
         ! ... only ionode reads the charge planes
         !
         IF ( ionode ) THEN
#if defined __HDF5
            CALL  read_rho_hdf5(h5desc , k,rho_plane)
#else
            CALL iotk_scan_dat( rhounit, "z" // iotk_index( k ), rho_plane )
#endif
         ENDIF
         !
         ! ... planes are sent to the destination processor (all processors in this image)
         !
         CALL mp_bcast( rho_plane, ionode_id, intra_image_comm )
         !
         IF( kowner(k) == me_group3 ) THEN
            !
            kk = k - fft_desc%my_i0r3p
            DO jj = 1, fft_desc%my_nr2p
               j = jj + fft_desc%my_i0r2p
               DO i = 1, nr1
                  rho(i+(jj-1)*nr1x+(kk-1)*ldr) = rho_plane(i+(j-1)*nr1)
               END DO
            END DO
            !
         END IF
         !
      END DO
      !
      DEALLOCATE( rho_plane )
      DEALLOCATE( kowner )
      !
      IF ( ionode ) THEN
         !
#if defined __HDF5
         CALL h5fclose_f(h5desc%file_id,ierr)
         DEALLOCATE ( h5desc)
#else
   !
         CALL iotk_scan_end( rhounit, "CHARGE-DENSITY" )
         !
         CALL iotk_close_read( rhounit )
#endif    
      END IF
      !
      RETURN
      !
    END SUBROUTINE read_rho_xml
    !
    !------------------------------------------------------------------------
    ! ... methods to write and read wavefunctions
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_wfc( iuni, ik, nk, kunit, ispin, nspin, wf0, ngw,   &
                          gamma_only, nbnd, igl, ngwl, filename, scalef, &
                          ionode, root_in_group, intra_group_comm,       &
                          inter_group_comm, parent_group_comm )
      !------------------------------------------------------------------------
      !
      USE mp_wave,    ONLY : mergewf
      USE control_flags,     ONLY : lwfnscf, lwfpbe0nscf  ! Lingzhu Kong
#if defined  __HDF5
      !USE hdf5_qe,    ONLY : evc_hdf5, read_data_hdf5, write_data_hdf5, &
      !                        evc_hdf5_write,  &
      !                       setup_file_property_hdf5, &
      !                       write_final_data, prepare_for_writing_final, &
      USE hdf5_qe                
      USE mp_global,    ONLY : inter_pool_comm, world_comm
      USE HDF5
#endif
      !
      IMPLICIT NONE
      !
      INTEGER,            INTENT(IN) :: iuni
      INTEGER,            INTENT(IN) :: ik, nk, kunit, ispin, nspin
      COMPLEX(DP),        INTENT(IN) :: wf0(:,:)
      INTEGER,            INTENT(IN) :: ngw
      LOGICAL,            INTENT(IN) :: gamma_only
      INTEGER,            INTENT(IN) :: nbnd
      INTEGER,            INTENT(IN) :: ngwl
      INTEGER,            INTENT(IN) :: igl(:)
      CHARACTER(LEN=256), INTENT(IN) :: filename
      REAL(DP),           INTENT(IN) :: scalef    
        ! scale factor, usually 1.0 for pw and 1/SQRT( omega ) for CP
      LOGICAL,            INTENT(IN) :: ionode
      INTEGER,            INTENT(IN) :: root_in_group, intra_group_comm, inter_group_comm, parent_group_comm
      !
      INTEGER                  :: j
      INTEGER                  :: iks, ike, ikt, igwx
      INTEGER                  :: ierr
      INTEGER                  :: ngroup, ipsour
      INTEGER,     ALLOCATABLE :: ipmask(:)
      INTEGER                  :: me_in_group, nproc_in_group, io_in_parent, nproc_in_parent, me_in_parent, my_group, io_group
#if defined __HDF5
      CHARACTER(LEN=256) :: filename_hdf5
#endif
      COMPLEX(DP), ALLOCATABLE :: wtmp(:)
      !
      ngroup          = mp_size( inter_group_comm )
      my_group        = mp_rank( inter_group_comm )
      me_in_group     = mp_rank( intra_group_comm )
      nproc_in_group  = mp_size( intra_group_comm )
      me_in_parent    = mp_rank( parent_group_comm )
      nproc_in_parent = mp_size( parent_group_comm )
      !
      ALLOCATE( ipmask( nproc_in_parent ) )
      !
      ! find out the group containing the ionode
      !
      io_group = 0
      IF( ionode ) io_group = my_group
      CALL mp_sum( io_group, parent_group_comm )
      !
      io_in_parent = 0
      IF( ionode ) io_in_parent = me_in_parent
      CALL mp_sum( io_in_parent, parent_group_comm )
      !
      CALL set_kpoints_vars( ik, nk, kunit, ngwl, igl, &
                             ngroup, ikt, iks, ike, igwx, ipmask, ipsour, &
                             ionode, root_in_group, intra_group_comm, inter_group_comm, parent_group_comm )
      !
      IF ( ionode ) THEN
#if defined  __HDF5
      filename_hdf5=trim(tmp_dir) //"evc.hdf5"
      CALL prepare_for_writing_final(evc_hdf5_write,inter_pool_comm,filename_hdf5,ik)
      CALL add_attributes_hdf5(evc_hdf5_write,ngw,"ngw",ik)
      CALL add_attributes_hdf5(evc_hdf5_write,gamma_only,"gamma_only",ik)
      CALL add_attributes_hdf5(evc_hdf5_write,igwx,"igwx",ik)
      CALL add_attributes_hdf5(evc_hdf5_write,nbnd,"nbnd",ik)
      CALL add_attributes_hdf5(evc_hdf5_write,ik,"ik",ik)
      CALL add_attributes_hdf5(evc_hdf5_write,nk,"nk",ik)
      CALL add_attributes_hdf5(evc_hdf5_write,ispin,"ispin",ik)
      CALL add_attributes_hdf5(evc_hdf5_write,nspin,"nspin",ik)
      CALL add_attributes_hdf5(evc_hdf5_write,scalef,"scale_factor",ik)
         !
#else
         CALL iotk_open_write( iuni, FILE = TRIM( filename ), ROOT="WFC", BINARY = .TRUE. )
         !
         CALL iotk_write_attr( attr, "ngw",          ngw, FIRST = .TRUE. )
         CALL iotk_write_attr( attr, "igwx",         igwx )
         CALL iotk_write_attr( attr, "gamma_only",   gamma_only )
         CALL iotk_write_attr( attr, "nbnd",         nbnd )
         CALL iotk_write_attr( attr, "ik",           ik )
         CALL iotk_write_attr( attr, "nk",           nk )
         CALL iotk_write_attr( attr, "ispin",        ispin )
         CALL iotk_write_attr( attr, "nspin",        nspin )
         CALL iotk_write_attr( attr, "scale_factor", scalef )
         !
         CALL iotk_write_empty( iuni, "INFO", attr )
#endif

         !
      END IF
      !
      ALLOCATE( wtmp( MAX( igwx, 1 ) ) )
      !
      wtmp = 0.0_DP
      ! Next 3 lines: Lingzhu Kong
      IF ( ( index(filename,'evc0') > 0 ) .and. (lwfnscf .or. lwfpbe0nscf) )THEN
         IF ( ionode ) OPEN(60,file='cp_wf.dat',status='unknown',form='unformatted')
      ENDIF

      DO j = 1, nbnd
         !
         IF ( ngroup > 1 ) THEN
            !
            IF( nk > 1 ) THEN
               IF ( ikt >= iks .AND. ikt <= ike ) &      
                  CALL mergewf( wf0(:,j), wtmp, ngwl, igl, me_in_group, &
                                nproc_in_group, root_in_group, intra_group_comm )
               !
               IF ( ipsour /= io_in_parent ) &
                  CALL mp_get( wtmp, wtmp, me_in_parent, &
                               io_in_parent, ipsour, j, parent_group_comm )
               !
            ELSE IF( my_group == io_group ) THEN
               !
               CALL mergewf( wf0(:,j), wtmp, ngwl, igl, &
                             me_in_group, nproc_in_group, root_in_group, intra_group_comm )
            END IF
            !
         ELSE
            !
            CALL mergewf( wf0(:,j), wtmp, ngwl, igl, &
                          me_in_parent, nproc_in_parent, io_in_parent, parent_group_comm )
            !
         END IF
         !
         IF ( ionode ) THEN

#if defined  __HDF5
            CALL write_evc(evc_hdf5_write,j,wtmp(1:igwx), ik)
#else
            CALL iotk_write_dat( iuni, "evc" // iotk_index( j ), wtmp(1:igwx) )
#endif
         ENDIF
         ! Next 3 lines : Lingzhu Kong
         IF ( ( index(filename,'evc0') > 0 ) .and. (lwfnscf .or. lwfpbe0nscf) ) THEN
            IF ( ionode ) write(60)wtmp(1:igwx) 
         ENDIF
         !
      END DO
      ! Next 4 lines : Lingzhu Kong
      IF ( ( index(filename,'evc0') > 0 ) .and. (lwfnscf .or. lwfpbe0nscf) )THEN
          IF ( ionode ) close(60)   !Lingzhu Kong
          write(*,*)'done writing evc0'
      ENDIF
      IF ( ionode ) then
#if defined __HDF5
         CALL h5fclose_f(evc_hdf5_write%file_id, ierr)
#else
         CALL iotk_close_write( iuni )
#endif
      endif
      !
      DEALLOCATE( wtmp )
      DEALLOCATE( ipmask )
      !
      RETURN
      !
    END SUBROUTINE write_wfc
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_wfc( iuni, ik, nk, kunit, ispin, &
                         nspin, wf, ngw, nbnd, igl, ngwl, filename, scalef, &
                         ionode, root_in_group, intra_group_comm, inter_group_comm, parent_group_comm, &
                         flink )
      !------------------------------------------------------------------------
      !
      USE mp_wave,   ONLY : splitwf

#if defined  __HDF5
      USE hdf5_qe
#endif
      !
      IMPLICIT NONE
      !
      INTEGER,            INTENT(IN)    :: iuni
      COMPLEX(DP),        INTENT(OUT)   :: wf(:,:)
      INTEGER,            INTENT(IN)    :: ik, nk
      INTEGER,            INTENT(IN)    :: kunit
      INTEGER,            INTENT(INOUT) :: ngw, nbnd, ispin, nspin
      INTEGER,            INTENT(IN)    :: ngwl
      INTEGER,            INTENT(IN)    :: igl(:)
      CHARACTER(LEN=256), INTENT(IN)    :: filename
      REAL(DP),           INTENT(OUT)   :: scalef
      LOGICAL,            INTENT(IN)    :: ionode
      INTEGER,            INTENT(IN)    :: root_in_group, intra_group_comm, inter_group_comm, parent_group_comm
      LOGICAL, OPTIONAL,  INTENT(IN)    :: flink
      !
      CHARACTER(LEN=256) :: filename_hdf5
      INTEGER                  :: j
      COMPLEX(DP), ALLOCATABLE :: wtmp(:)
      INTEGER                  :: ierr
      INTEGER                  :: iks, ike, ikt
      INTEGER                  :: igwx, igwx_, ik_, nk_
      INTEGER                  :: ngroup, ipdest
      INTEGER,     ALLOCATABLE :: ipmask(:)
      LOGICAL                  :: flink_
      INTEGER                  :: me_in_group, nproc_in_group, io_in_parent, nproc_in_parent, me_in_parent, my_group, io_group
      !
      flink_ = .FALSE.
      IF( PRESENT( flink ) ) flink_ = flink
      !
      ngroup          = mp_size( inter_group_comm )
      my_group        = mp_rank( inter_group_comm )
      me_in_group     = mp_rank( intra_group_comm )
      nproc_in_group  = mp_size( intra_group_comm )
      me_in_parent    = mp_rank( parent_group_comm )
      nproc_in_parent = mp_size( parent_group_comm )
      !
      ALLOCATE( ipmask( nproc_in_parent ) )
      !
      ! find out the group containing the ionode
      !
      io_group = 0
      IF( ionode ) io_group = my_group
      CALL mp_sum( io_group, parent_group_comm )
      !
      ! find out the io task id in the parent group
      !
      io_in_parent = 0
      IF( ionode ) io_in_parent = me_in_parent
      CALL mp_sum( io_in_parent, parent_group_comm )
      !
      CALL set_kpoints_vars( ik, nk, kunit, ngwl, igl, &
                             ngroup, ikt, iks, ike, igwx, ipmask, ipdest, &
                             ionode, root_in_group, intra_group_comm, inter_group_comm, parent_group_comm )
      !
      !  if flink = .true. we are following a link and the file is
      !  already opened for read
      !
      ierr = 0
      !
#if !defined __HDF5
      IF ( ionode .AND. .NOT. flink_ ) &
         CALL iotk_open_read( iuni, FILE = filename, &
                              BINARY = .TRUE., IERR = ierr )
      !
      CALL mp_bcast( ierr, io_in_parent, parent_group_comm )
      !
      CALL errore( 'read_wfc ', &
                   'cannot open restart file for reading', ierr )
#endif
      !
      IF ( ionode ) THEN
          !
#if defined  __HDF5
          !filename_hdf5=trim(tmp_dir) //"evc.hdf5"
          filename_hdf5=filename
          CALL prepare_for_reading_final(evc_hdf5_write,evc_hdf5_write%comm,filename_hdf5,ik)
          CALL read_attributes_hdf5(evc_hdf5_write,ngw,"ngw",ik)
          CALL read_attributes_hdf5(evc_hdf5_write,nbnd,"nbnd",ik)
          CALL read_attributes_hdf5(evc_hdf5_write,ik_,"ik",ik)
          CALL read_attributes_hdf5(evc_hdf5_write,nk_,"ik",ik)
          CALL read_attributes_hdf5(evc_hdf5_write,ispin,"ispin",ik)
          CALL read_attributes_hdf5(evc_hdf5_write,nspin,"nspin",ik)
          CALL read_attributes_hdf5(evc_hdf5_write,igwx_,"igwx",ik)
          CALL read_attributes_hdf5(evc_hdf5_write,scalef,"scale_factor",ik)
#else
          CALL iotk_scan_empty( iuni, "INFO", attr )
          !
          CALL iotk_scan_attr( attr, "ngw",          ngw )
          CALL iotk_scan_attr( attr, "nbnd",         nbnd )
          CALL iotk_scan_attr( attr, "ik",           ik_ )
          CALL iotk_scan_attr( attr, "nk",           nk_ )
          CALL iotk_scan_attr( attr, "ispin",        ispin )
          CALL iotk_scan_attr( attr, "nspin",        nspin )
          CALL iotk_scan_attr( attr, "igwx",         igwx_ )
          CALL iotk_scan_attr( attr, "scale_factor", scalef )
          !
#endif

      END IF
      !
      CALL mp_bcast( ngw,    io_in_parent, parent_group_comm )
      CALL mp_bcast( nbnd,   io_in_parent, parent_group_comm )
      CALL mp_bcast( ik_,    io_in_parent, parent_group_comm )
      CALL mp_bcast( nk_,    io_in_parent, parent_group_comm )
      CALL mp_bcast( ispin,  io_in_parent, parent_group_comm )
      CALL mp_bcast( nspin,  io_in_parent, parent_group_comm )
      CALL mp_bcast( igwx_,  io_in_parent, parent_group_comm )
      CALL mp_bcast( scalef, io_in_parent, parent_group_comm )
      !
      ALLOCATE( wtmp( MAX( igwx_, igwx ) ) )
      !
      DO j = 1, nbnd
         !
         IF ( j <= SIZE( wf, 2 ) ) THEN
            !
            IF ( ionode ) THEN 
               !
#if defined __HDF5
             CALL read_evc(evc_hdf5_write,j,wtmp(1:igwx_),ik)
             !  CALL iotk_scan_dat( iuni, &
             !                      "evc" // iotk_index( j ), wtmp(1:igwx_) )
#else
               CALL iotk_scan_dat( iuni, &
                                   "evc" // iotk_index( j ), wtmp(1:igwx_) )
 
#endif
               !
               IF ( igwx > igwx_ ) wtmp((igwx_+1):igwx) = 0.0_DP
               ! ===========================================================
               !       Lingzhu Kong
               !IF ( j .eq. 1)write(*,'(10f12.5)')(wtmp(i),i=1,igwx_)
               ! ===========================================================
               !
            END IF
            !
            IF ( ngroup > 1 ) THEN
               !
               IF( nk_ > 1 ) THEN
                  !
                  IF ( ipdest /= io_in_parent ) &
                     CALL mp_put( wtmp, wtmp, me_in_parent, &
                               io_in_parent, ipdest, j, parent_group_comm )
                  !
                  IF ( ( ikt >= iks ) .AND. ( ikt <= ike ) ) &
                     CALL splitwf( wf(:,j), wtmp, ngwl, igl, me_in_group, &
                                nproc_in_group, root_in_group, intra_group_comm )
                  !
               ELSE IF( my_group == io_group ) THEN

                  CALL splitwf( wf(:,j), wtmp, ngwl, igl, &
                             me_in_group, nproc_in_group, root_in_group, intra_group_comm )
               END IF
               !
            ELSE
               !
               CALL splitwf( wf(:,j), wtmp, ngwl, igl, &
                             me_in_parent, nproc_in_parent, io_in_parent, parent_group_comm )
               !
            END IF
            !
         END IF
         !
      END DO
      !
#if !defined __HDF5
      IF ( ionode .AND. .NOT. flink_ ) CALL iotk_close_read( iuni )
#endif
      !
      DEALLOCATE( wtmp )
      DEALLOCATE( ipmask )
      !
      RETURN
      !
    END SUBROUTINE read_wfc
    !        
END MODULE oldxml_xml_io_base
