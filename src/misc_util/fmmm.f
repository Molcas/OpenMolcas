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
      SUBROUTINE FMMM(A,B,C,NROW,NCOL,NSUM)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NROW,NSUM),B(NSUM,NCOL),C(NROW,NCOL)
      PARAMETER (KS=48)
C
      DO 100 I = 1, NROW
         DO 200 J = 1, NCOL
            C(I,J) = 0.0D0
 200     CONTINUE
 100  CONTINUE
C
      DO 300 KK = 1, NSUM, KS
         DO 400 I = 1, NROW
            DO 500 J = 1, NCOL
               T = C(I,J)
               DO 600 K = KK,MIN(KK+KS-1,NSUM)
                  T = T + B(K,J) * A(I,K)
 600           CONTINUE
               C(I,J) = T
 500        CONTINUE
 400     CONTINUE
 300  CONTINUE
C
      RETURN
      END
