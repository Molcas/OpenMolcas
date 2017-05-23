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
      SUBROUTINE DNAXPY(N,M,A,INCA,X,INCXI,INCXO,Y,INCYI,INCYO)
C
C     MULTIPLY A VECTOR, X, BY A SCALAR, ADD TO A VECTOR, Y, AND
C     STORE THE RESULT IN THE VECTOR Y. REPEAT THIS N TIMES.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A((N-1)*INCA+1),
     1       X(((M-1)*INCXI+1)*((N-1)*INCXO+1)),
     2       Y(((M-1)*INCYI+1)*((N-1)*INCYO+1))
      DO 10 I = 1, N
         CALL DAXPY_(M,A(1+(I-1)*INCA),
     1                X(1+(I-1)*INCXO),INCXI,
     2                Y(1+(I-1)*INCYO),INCYI)
   10 CONTINUE
      RETURN
      END
