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
      SUBROUTINE DYAX(N,ALPHA,X,INCX,Y,INCY)
C
C     MULTIPLY A VECTOR, X, BY A SCALAR AND STORE THE RESULT IN
C     THE VECTOR Y.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(1+(N-1)*INCX), Y(1+(N-1)*INCY)
      DO 10 I = 1, N
         Y(1+(I-1)*INCY) = ALPHA * X(1+(I-1)*INCX)
   10 CONTINUE
      RETURN
      END
