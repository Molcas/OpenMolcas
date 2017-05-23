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
      FUNCTION CHO_DSUMELM(VEC,N)
C
C     Purpose: sum elements of double precision vector.
C
      IMPLICIT NONE
      REAL*8  CHO_DSUMELM
      REAL*8  VEC(*), DSUM
      INTEGER I, N

      IF (N .GT. 0) THEN
         DSUM = VEC(1)
         DO I = 2,N
            DSUM = DSUM + VEC(I)
         END DO
      ELSE
         DSUM = 0.0D0
      END IF

      CHO_DSUMELM = DSUM

      END
