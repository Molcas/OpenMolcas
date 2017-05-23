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
      SUBROUTINE CHO_INVPCK(IJ,I,J,LOW)
C
C     Purpose: invert triangular packing index ij to
C              rectangular indices i and j.
C              Flag LOW specifies packing convention:
C              LOW = T: i>=j
C              LOW = F: i<=j
C
#include "implicit.fh"
      LOGICAL LOW

      PARAMETER (ONE = 1.0D0, TWO = 2.0D0, THREE = 3.0D0, EIGHT = 8.0D0)

      IF (IJ .GT. 0) THEN

         XX = EIGHT*DBLE(IJ) - THREE
         XI = (ONE + SQRT(XX))/TWO
         I  = INT(XI)
         J  = IJ - I*(I - 1)/2

         IF (.NOT. LOW) THEN
            ITMP = I
            I = J
            J = ITMP
         END IF

      ELSE

         I = -1
         J = -2

      END IF

      END
