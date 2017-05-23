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
      INTEGER FUNCTION CHO_ISUMELM(IVEC,N)
C
C     Purpose: sum elements of integer vector.
C
      IMPLICIT NONE
      INTEGER IVEC(*)
      INTEGER I, N, ISUM

      IF (N .GT. 0) THEN
         ISUM = IVEC(1)
         DO I = 2,N
            ISUM = ISUM + IVEC(I)
         END DO
      ELSE
         ISUM = 0
      END IF

      CHO_ISUMELM = ISUM

      END
