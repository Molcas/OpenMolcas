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
      INTEGER FUNCTION CHO_IDOT(N,IVEC1,INC1,IVEC2,INC2)
      IMPLICIT NONE
      INTEGER N, INC1, INC2
      INTEGER IVEC1(*), IVEC2(*)

      INTEGER I, ISUM, I1, I2

      CHO_IDOT = 0
      IF (N.LT.1) RETURN
      ISUM = 0
      IF (INC1.EQ.1 .AND. INC2.EQ.1) THEN
         DO I = 1,N
            ISUM = ISUM + IVEC1(I)*IVEC2(I)
         END DO
      ELSE
         I1 = 1
         I2 = 1
         IF (INC1 .LT. 0) I1 = (-N+1)*INC1 + 1
         IF (INC2 .LT. 0) I2 = (-N+1)*INC2 + 1
         DO I = 1,N
            ISUM = ISUM + IVEC1(I1)*IVEC2(I2)
            I1 = I1 + INC1
            I2 = I2 + INC2
         END DO
      END IF
      CHO_IDOT = ISUM

      END
