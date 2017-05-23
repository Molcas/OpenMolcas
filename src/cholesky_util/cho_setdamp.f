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
      SUBROUTINE CHO_SETDAMP()
C
C     Purpose: set screening damping, unless user-defined.
C
      IMPLICIT NONE
#include "cholesky.fh"
      INTEGER I

      DO I = 1,2
         IF (DAMP(I) .LT. 0.0D0) THEN
            IF (THRCOM .GT. 9.99D-3) THEN ! >= 1.0D-2
               DAMP(I) = 1.0D7
            ELSE IF (THRCOM .GT. 9.99D-4) THEN ! >= 1.0D-3
               DAMP(I) = 1.0D6
            ELSE IF (THRCOM .GT. 9.99D-5) THEN ! >= 1.0D-4
               DAMP(I) = 1.0D5
            ELSE IF (THRCOM .GT. 9.99D-6) THEN ! >= 1.0D-5
               DAMP(I) = 1.0D4
            ELSE IF (THRCOM .GT. 9.99D-7) THEN ! >= 1.0D-6
               DAMP(I) = 1.0D3
            ELSE IF (THRCOM .GT. 9.99D-8) THEN ! >= 1.0D-7
               DAMP(I) = 1.0D2
            ELSE IF (THRCOM .GT. 9.99D-9) THEN ! >= 1.0D-8
               DAMP(I) = 1.0D1
            ELSE
               DAMP(I) = 1.0D0
            END IF
         END IF
      END DO

      END
