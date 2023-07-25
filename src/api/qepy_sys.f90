! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE qepy_sys
   USE kinds,               ONLY : DP
   !
   IMPLICIT NONE
   PUBLIC
   !
   CHARACTER(LEN=512) :: command_line = ' '
   !
#if defined(__MPI)
   logical :: is_mpi = .TRUE.
#else
   logical :: is_mpi = .FALSE.
#endif
#if defined(_OPENMP)
   logical :: is_openmp = .TRUE.
#else
   logical :: is_openmp = .FALSE.
#endif
   !
   INTERFACE COMMAND_ARGUMENT_COUNT
      MODULE PROCEDURE qepy_my_iargc
   END INTERFACE
   !
   INTERFACE GET_COMMAND_ARGUMENT
      MODULE PROCEDURE qepy_my_getarg
   END INTERFACE
   !
CONTAINS
   !
   INTEGER FUNCTION qepy_my_iargc()
      IMPLICIT NONE

      CHARACTER(LEN=1) :: previous, current
      INTEGER :: i
      !IF (LEN_TRIM(command_line) == 0) THEN
         !qepy_my_iargc = COMMAND_ARGUMENT_COUNT()
         !RETURN
      !ENDIF

      qepy_my_iargc = 0
      previous = ' '
      DO i=1,LEN_TRIM(command_line)
         current = command_line(i:i)
         IF ( current /= ' ' .AND. previous == ' ' ) qepy_my_iargc = qepy_my_iargc+1
         previous = current
      END DO

   END FUNCTION qepy_my_iargc
   !
   SUBROUTINE qepy_my_getarg( narg, arg )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: narg 
      CHARACTER(LEN=*), INTENT(OUT) :: arg
      CHARACTER(LEN=1) :: previous, current
      INTEGER :: iarg, i, indx

      !IF (LEN_TRIM(command_line) == 0) THEN
         !CALL GET_COMMAND_ARGUMENT(narg, arg)
         !RETURN
      !ENDIF

      iarg = 0
      previous = ' '
      arg = ' '
      indx= 0
      DO i=1,LEN_TRIM(command_line)
         current = command_line(i:i)
         IF ( current /= ' ' .AND. previous == ' ' ) iarg = iarg+1
         IF ( iarg == narg ) THEN
            indx = indx + 1
            arg(indx:indx) = current           
            IF ( indx == LEN(arg) ) RETURN
         ELSE IF ( iarg > narg ) THEN
            RETURN
         END IF
         previous = current
      END DO

   END SUBROUTINE qepy_my_getarg 
END MODULE qepy_sys
