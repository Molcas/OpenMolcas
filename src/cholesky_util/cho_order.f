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
      SUBROUTINE CHO_ORDER(VEC,LVEC,IJOB)
C
C     Purpose: sort elements of VEC in
C              IJOB = -1: descending order
C              IJOB =  1: ascending  order
C              (all other IJOB values are ignored)
C
#include "implicit.fh"
      DIMENSION VEC(LVEC)

      IF (IJOB .EQ. -1) THEN

         DO I = 1,LVEC-1
            VMAX = VEC(I)
            IMAX = I
            DO J = I+1,LVEC
               IF (VEC(J) .GT. VMAX) THEN
                  VMAX = VEC(J)
                  IMAX = J
               END IF
            END DO
            IF (IMAX .NE. I) THEN
               VEC(IMAX) = VEC(I)
               VEC(I)    = VMAX
            END IF
         END DO

      ELSE IF (IJOB .EQ. 1) THEN

         DO I = 1,LVEC-1
            VMIN = VEC(I)
            IMIN = I
            DO J = I+1,LVEC
               IF (VEC(J) .LT. VMIN) THEN
                  VMIN = VEC(J)
                  IMIN = J
               END IF
            END DO
            IF (IMIN .NE. I) THEN
               VEC(IMIN) = VEC(I)
               VEC(I)    = VMIN
            END IF
         END DO

      END IF

      END
