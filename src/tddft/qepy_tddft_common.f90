!
! Copyright (C) 2001-2019 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
MODULE qepy_tddft_common
  !----------------------------------------------------------------------
  !  ... Compute optical absorption spectrum by real-time TDDFT 
  !  ... References:
  !      (1) Phys. Rev. B 73, 035408 (2006)
  !      (2) http://www.netlib.org/linalg/html_templates/Templates.html
  !                                             Xiaofeng Qian, MIT (2008)
  !----------------------------------------------------------------------
  USE kinds,                       ONLY : dp
  !
  IMPLICIT NONE

  !-- tddft variables ----------------------------------------------------
  complex(dp), allocatable :: tddft_psi(:,:,:), b(:,:)
  complex(dp), allocatable :: tddft_hpsi(:,:), tddft_spsi(:,:)
  complex(dp), allocatable :: tddft_Ppsi(:,:)        ! PAW correction to forces (Ehrenfest)
  real(dp), allocatable :: charge(:), dipole(:,:), quadrupole(:,:,:)
  complex(dp), allocatable :: circular(:,:), circular_local(:)

END MODULE qepy_tddft_common
