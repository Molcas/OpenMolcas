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
      SUBROUTINE CHO_DZERO(VEC,N)
C
C     Purpose: zero double precision vector.
C
      IMPLICIT NONE
      REAL*8  VEC(*)
      INTEGER I, N

      DO I = 1,N
         VEC(I) = 0.0D0
      END DO

      END
