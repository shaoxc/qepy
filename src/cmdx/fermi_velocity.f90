!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
!
! Usage :
! $ fermi_velocity.x -in {pw.x input file}
! Then it generates vfermi.frmsf (for nspin = 1, 4) or
! vfermi1.frmsf and vfermi2.frmsf (for nspin = 2)
!
!----------------------------------------------------------------------------
SUBROUTINE fermi_velocity()
  !--------------------------------------------------------------------------
  !
  USE input_parameters,     ONLY : prefix, outdir
  USE io_files,             ONLY : prefix_ => prefix, tmp_dir
  USE mp_global,            ONLY : mp_startup
  USE environment,          ONLY : environment_start, environment_end
  USE read_input,           ONLY : read_input_file
  USE command_line_options, ONLY : input_file_
  USE kinds,                ONLY : DP
  USE wvfct,                ONLY : nbnd, et
  USE start_k,              ONLY : nk1, nk2, nk3
  USE cell_base,            ONLY : at, alat
  USE lsda_mod,             ONLY : nspin
  USE ener,                 ONLY : ef, ef_up, ef_dw
  USE klist,                ONLY : nks, two_fermi_energies
  USE fermisurfer_common,   ONLY : b_low, b_high, rotate_k_fs, write_fermisurfer
  USE constants,            ONLY : tpi
  !
  IMPLICIT NONE
  !
  INTEGER :: i1, i2, i3, ibnd, ii, ikp(3), ikm(3), ispin, ns, nk
  REAL(DP) :: de(3), ef1, ef2
  INTEGER,ALLOCATABLE :: equiv(:,:,:)
  REAL(DP),ALLOCATABLE :: eig(:,:,:,:,:), vf(:,:,:,:,:)
  LOGICAL :: needwf = .FALSE.
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  CALL mp_startup ()
  CALL environment_start ('FERMI_VELOCITY')
  !
  ! ... Read pw.x input file and get prefix and outdir
  !
  CALL read_input_file ('PW', input_file_)
  !
  prefix_ = TRIM(prefix)
  tmp_dir = trimcheck(outdir)
  !
  ! ... Read XML file generated by pw.x
  !
  CALL read_file_new( needwf)
  !
  ! ... Number of k and spin for each magnetic treatment
  !
  IF (nspin == 2) THEN
     ns = 2
     IF(two_fermi_energies) THEN
        ef1 = ef_up
        ef2 = ef_dw
     ELSE
        ef1 = ef
        ef2 = ef
     END IF
  ELSE
     ns = 1
  END IF
  nk = nks / ns
  !
  ! ... Find equivalent k point in irr-BZ for whole BZ
  !
  ALLOCATE(equiv(nk1, nk2, nk3))
  CALL rotate_k_fs(equiv)
  !
  ALLOCATE(vf(b_low:b_high, nk1, nk2, nk3, ns), &
  &       eig(b_low:b_high, nk1, nk2, nk3, ns))
  !
  ! ... Map e_k into whole BZ (Measured from E_F)
  !
  DO i3 = 1, nk3
     DO i2 = 1, nk2
        DO i1 = 1, nk1
           IF(nspin == 2) THEN
              eig(b_low:b_high,i1,i2,i3,1) = et(b_low:b_high, equiv(i1,i2,i3)     ) - ef1
              eig(b_low:b_high,i1,i2,i3,2) = et(b_low:b_high, equiv(i1,i2,i3) + nk) - ef2
           ELSE
              eig(b_low:b_high,i1,i2,i3,1) = et(b_low:b_high, equiv(i1,i2,i3)     ) - ef
           END IF
        END DO
     END DO
  END DO
  !
  ! ... Compute Fermi velocity in the atomic unit
  !
  DO i3 = 1, nk3
     DO i2 = 1, nk2
        DO i1 = 1, nk1
           !
           DO ispin = 1, ns
              !
              DO ibnd = b_low, b_high
                 !
                 DO ii = 1, 3
                    !
                    ikp(1:3) = (/i1, i2, i3/) - 1
                    ikp(ii) = ikp(ii) + 1
                    ikp(1:3) = MODULO(ikp(1:3), (/nk1, nk2, nk3/)) + 1
                    !
                    ikm(1:3) = (/i1, i2, i3/) - 1
                    ikm(ii) = ikm(ii) - 1
                    ikm(1:3) = MODULO(ikm(1:3), (/nk1, nk2, nk3/)) + 1
                    !
                    de(ii) = eig(ibnd, ikp(1), ikp(2), ikp(3), ispin) &
                    &      - eig(ibnd, ikm(1), ikm(2), ikm(3), ispin)
                    !
                 END DO
                 !
                 de(1:3) = 0.5_dp * de(1:3) * REAL((/nk1, nk2, nk3/), DP)
                 de(1:3) = matmul(at(1:3,1:3), de(1:3)) * alat / tpi
                 vf(ibnd, i1, i2, i3, ispin) = SQRT(DOT_PRODUCT(de, de))
                 !
              END DO ! ibnd = 1, nbnd
              !
           END DO ! ispin = 1, ns
        END DO ! i1 = 1, nk1
     END DO ! i2 = 1, nk2
  END DO ! i3 = 1, nk3
  !
  ! ... Output in the FermiSurfer format
  !
  IF (nspin == 2) THEN
     CALL write_fermisurfer(eig(b_low:b_high, 1:nk1, 1:nk2, 1:nk3, 1), &
     &                       vf(b_low:b_high, 1:nk1, 1:nk2, 1:nk3, 1), "vfermi1.frmsf")
     CALL write_fermisurfer(eig(b_low:b_high, 1:nk1, 1:nk2, 1:nk3, 2), &
     &                       vf(b_low:b_high, 1:nk1, 1:nk2, 1:nk3, 2), "vfermi2.frmsf")
  ELSE
     CALL write_fermisurfer(eig(b_low:b_high, 1:nk1, 1:nk2, 1:nk3, 1), &
     &                       vf(b_low:b_high, 1:nk1, 1:nk2, 1:nk3, 1), "vfermi.frmsf")
  END IF
  !
  DEALLOCATE(vf, eig, equiv)
  !
  CALL environment_end ('FERMI_VELOCITY')
  CALL stop_pp
  !
CONTAINS
!


END SUBROUTINE fermi_velocity