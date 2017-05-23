************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE CHO_GETMAXSHL(DIASH,SMAX,ISHLAB)
C
C     Purpose: Get max. shell pair and update DIASH.
C
      IMPLICIT NONE
      REAL*8  DIASH(*)
      REAL*8  SMAX
      INTEGER ISHLAB
#include "cholesky.fh"

      INTEGER JSHLAB

      CHARACTER*13 SECNAM
      PARAMETER (SECNAM = 'CHO_GETMAXSHL')

      SMAX   = -1.0D9
      ISHLAB = -1
      DO JSHLAB = 1,NNSHL
         IF (DIASH(JSHLAB) .GT. SMAX) THEN
            SMAX   = DIASH(JSHLAB)
            ISHLAB = JSHLAB
         END IF
      END DO

      IF (ISHLAB .LT. 1) THEN
         CALL CHO_QUIT('Error in '//SECNAM,104)
      ELSE
         DIASH(ISHLAB) = 0.0D0
      END IF

      END
