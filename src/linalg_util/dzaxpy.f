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
      SUBROUTINE DZAXPY(N,DA,DX,INCX,DY,INCY,DZ,INCZ)
C
C     MULTIPLY A VECTOR, X, BY A SCALAR, ADD TO A VECTOR, Y, AND
C     STORE THE RESULT IN THE VECTOR Z.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DX(*),DY(*),DZ(*)
      DATA ZERO /0.0D+00/
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IY = 1
      IZ = 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      IF(INCZ.LT.0)IZ = (-N+1)*INCZ + 1
      IF (DA .NE. ZERO ) THEN
         IX = 1
         IF(INCX.LT.0)IX = (-N+1)*INCX + 1
         DO I = 1,N
           DZ(IY) = DY(IY) + DA*DX(IX)
           IX = IX + INCX
           IY = IY + INCY
           IZ = IZ + INCZ
         End Do
      ELSE
         DO I = 1,N
           DZ(IY) = DY(IY)
           IY = IY + INCY
           IZ = IZ + INCZ
         End Do
      END IF
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,4)
      IF (DA .NE. ZERO ) THEN
         IF( M .NE. 0 ) THEN
            DO I = 1,M
              DZ(I) = DY(I) + DA*DX(I)
            END DO
            IF( N .LT. 4 ) RETURN
         END IF
         MP1 = M + 1
         DO I = MP1,N,4
           DZ(I) = DY(I) + DA*DX(I)
           DZ(I + 1) = DY(I + 1) + DA*DX(I + 1)
           DZ(I + 2) = DY(I + 2) + DA*DX(I + 2)
           DZ(I + 3) = DY(I + 3) + DA*DX(I + 3)
         END DO
      ELSE
         IF( M .NE. 0 ) THEN
            DO I = 1,M
              DZ(I) = DY(I)
            END DO
            IF( N .LT. 4 ) RETURN
         END IF
         MP1 = M + 1
         DO I = MP1,N,4
           DZ(I) = DY(I)
           DZ(I + 1) = DY(I + 1)
           DZ(I + 2) = DY(I + 2)
           DZ(I + 3) = DY(I + 3)
         END DO
      END IF
      RETURN
      END
