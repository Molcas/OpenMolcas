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
      INTEGER FUNCTION IDMIN (N, X, INCX)
C
C     FINDS THE INDEX OF ELEMENT HAVING MIN. VALUE.
C
      INTEGER    N, INCX
      REAL*8     X(*)
      INTEGER    I, IXX
      REAL*8     MIN
C
      IDMIN = 0
      IF (N .GE. 1) THEN
         IDMIN = 1
         IF (N .NE. 1) THEN
            IF (INCX .NE. 1) THEN
C                                  CODE FOR INCREMENT NOT EQUAL TO 1
               IXX = 1
               MIN = X(1)
               IXX = IXX + INCX
               DO 10  I=2, N
                  IF (X(IXX) .LT. MIN) THEN
                     IDMIN = I
                     MIN = X(IXX)
                  END IF
                  IXX = IXX + INCX
   10          CONTINUE
            ELSE
C                                  CODE FOR INCREMENT EQUAL TO 1
               MIN = X(1)
               DO 20  I=2, N
                  IF (X(I) .LT. MIN) THEN
                     IDMIN = I
                     MIN = X(I)
                  END IF
   20          CONTINUE
            END IF
         END IF
      END IF
      RETURN
      END
