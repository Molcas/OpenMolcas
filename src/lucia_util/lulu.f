************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1988, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE LULU(A,L,U,NDIM)
C
C LU DECOMPOSITION OF MATRIX A
C
C     A = L * U
C
C WHERE L IS A LOWER TRIANGULAR MATRIX WITH A
C UNIT DIAGONAL AND U IS AN UPPER DIAGONAL
C
C L AND U ARE STORED AS ONE DIMENSIONAL ARRAYS
C
C   L(I,J) = L(I*(I-1)/2 + J ) ( I .GE. J )
C
C   U(I,J) = U(J*(J-1)/2 + I ) ( J .GE. I )
C
C THIS ADRESSING SCHEMES SUPPORTS VECTORIZATION OVER COLUMNS
C FOR L AND  OVER ROWS FOR U .
C
C
C NO PIVOTING IS DONE HERE , SO THE SCHEME GOES :
C
C     LOOP OVER R=1, NDIM
C        LOOP OVER J = R, NDIM
C          U(R,J) = A(R,J) - SUM(K=1,R-1) L(R,K) * U(K,J)
C        END OF LOOP OVER J
C
C        LOOP OVER I = R+1, NDIM
C          L(I,R) = (A(I,R) - SUM(K=1,R-1)L(I,K) * U(K,R) ) /U(R,R)
C        END OF LOOP OVER I
C     END OF LOOP OVER R
C
C JEPPE OLSEN , OCTOBER 1988
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NDIM,NDIM)
      REAL*8  L(*),U(*)
      REAL * 8  INPROD
      INTEGER R
C
C
      DO 1000 R = 1, NDIM
C
        DO 100 J = R, NDIM
         U(J*(J-1)/2 + R ) = A(R,J) -
     &   INPROD(L(R*(R-1)/2+1),U(J*(J-1)/2+1),R-1)
  100   CONTINUE
C
        XFACI = 1.0D0/ U(R*(R+1)/2)
        L(R*(R+1)/2 ) = 1.0D0
        DO 200 I = R+1, NDIM
          L(I*(I-1)/2 + R) = (A(I,R) -
     &   INPROD(L(I*(I-1)/2+1),U(R*(R-1)/2+1),R-1) ) * XFACI
  200  CONTINUE
C
 1000 CONTINUE
C
      NTEST = 0
      IF ( NTEST .NE. 0 ) THEN
         WRITE(6,*) ' L MATRIX '
         CALL PRSYM(L,NDIM)
         WRITE(6,*) ' U MATRIX ( TRANSPOSED ) '
         CALL PRSYM(U,NDIM)
      END IF
C
      RETURN
      END
*
