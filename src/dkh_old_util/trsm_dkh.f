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
      SUBROUTINE TRSM_DKH(A,B,C,N,H,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N*(N+1)/2),B(N,N),C(N*(N+1)/2),H(N,N),W(N,N)
C
C     TRANSFORM SYMMETRIC MATRIX A BY UNITARY TRANSFORMATION
C     IN B. RESULT IS IN C
C
      IJ=0
      DO 10 I=1,N
         DO 11 J=1,I
            IJ=IJ+1
            C(IJ)=0.D0
            W(I,J)=A(IJ)
            W(J,I)=A(IJ)
            H(I,J)=0.D0
            H(J,I)=0.D0
11       CONTINUE
10    CONTINUE
      DO 2 I=1,N
      DO 3 L=1,N
      DO 4 K=1,N
      H(I,L)=B(K,I)*W(K,L)+H(I,L)
4     CONTINUE
3     CONTINUE
2     CONTINUE
      IJ=0
      DO 21 I=1,N
      DO 22 J=1,I
      IJ=IJ+1
      DO 23 L=1,N
      C(IJ)=H(I,L)*B(L,J)+C(IJ)
23    CONTINUE
22    CONTINUE
21    CONTINUE
      RETURN
      END
