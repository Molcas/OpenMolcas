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
* Copyright (C) 1988, H. Rieger                                        *
************************************************************************
C**********************************************************************
C
C        PROGRAMMBIBLIOTHEK RHRZ BONN        23/11/88       DCOPIV
C                                            FORTRAN 77     IBM 3081
C
C
C NAME:    DCOPIV
C
C PURPOSE:
C
C This program calculates the solution of the system of linear
C equations  A*X = B  for several right-hand sides. The method used is
C complete pivoting which is generally considered a stable method,
C even if the matrix A is ill-conditioned.
C A is decomposed into the product of a lower triangular matrx L
C with main diagonal elements equal to 1 and an upper triangular
C matrix U such that:  L*U = A. While L is not saved, U is over-
C written on A.
C If the matrix A is nonsingular, the solutions are stored in B.
C
C USAGE:   CALL DCOPIV(A,B,N,M,D,EPS,DET,EX,CTR,S)
C
C PARAMETERS:
C
C N,M:     INTEGER, are the order of the matrix A (N) and the
C          number of right-hand sides (M).
C
C D:       INTEGER, is the actual dimension of the arrays as they
C          are declared in the main program. D >= N is necessary.
C
C A:       REAL*8, declared as DIMENSION A(1:D,1:N), is the coefficient
C          matrix. On return it will be over-written by U.
C
C B:       REAL*8, declared as DIMENSION B(1:D,1:M), contains the right
C          hand sides (columnwise). On return each column occupies the
C          corresponding solution.
C
C EPS:     REAL*8, is the tolerance: if a pivotal element is
C          (absolutely) smaller than EPS, matrix A is considered
C          to be singular.
C
C DET:     REAL*8, will contain the determinant (DET = 0 is set,
C          if A is singular). It is scaled such that
C          1.0D-10 <= ABS(DET) <= 1.0D10 .
C
C EX:      INTEGER, is the scale factor for DET. The determinant
C          of A can be obtained from  DET*10**EX .
C
C CTR:     INTEGER, is a control variable with the following
C          meaning on input:
C          CTR < 0: evaluation of just the determinant. Parameter
C                   B is not used at all.
C          CTR = 0: solution of the system of equations and evaluation
C                   of the determinant.
C          CTR > 0: solution of the system of equations only. EX = 0
C                   and DET = 0 (1) will be set, if A is singular
C                   (regular).
C          On return it contains the error code:
C          CTR = -1: parameter fault (see REMARK 1).
C          CTR =  0: program worked.
C          CTR =  1: matrix A is singular.
C
C S:       INTEGER, declared as DIMENSION S(1:N), is used for
C          auxiliary storage only and has no meaning otherwise.
C
C REMARKS: (1) If N < 1 or M < 1, CTR = -1 is set and control is
C              returned to the calling program.
C          (2) The original contents of matrices A and B are
C              destroyed on return.
C          (3) Choosing B = I (unit matrix) with M = N will
C              generate the inverse matrix of A in B.
C          (4) This program makes no use of any special forms or
C              properties of the coefficient matrix (Hessenberg or
C              tridiagonal form, band matrices, symmetric or sparse
C              matrices). Special algorithms should be used in
C              that case, if the order is sufficiently high.
C
C REF.:    J.H. Wilkinson,
C          The Algebraic Eigenvalue Problem,
C          Oxford University Press, pp 212-214, 1965.
C
C AUTHOR:        H. RIEGER, RHRZ
C INSTALLATION:  IBM 3081, VS FORTRAN V2 (OPT=3)
C
C ACCESS:
C                --- MVS ---              --- VM/CMS ---
C LOAD MODULE:   SYS3.FORTLIB(DCOPIV)     FORTLIB  TXTLIB P (DCOPIV)
C DESCRIPTION:   SYS3.INFOLIB(DCOPIV)     INFOLIB  MACLIB P (DCOPIV)
C SOURCE MODULE: SYS3.SYMLIB.FORTRAN(DCOPIV)
C
C**********************************************************************
C
      SUBROUTINE   DCOPIV(A,B,N,M,iD,EPS,DET,iEX,iCTR,iS)
      implicit real*8(a-h,o-z)
      DIMENSION    A(iD,N),B(iD,M),iS(N)
C
      IF (N .LT. 1 .OR. M .LT. 1) GO TO 7
      iEX = 0
      DET = 1.0D0
      IF (N .EQ. 1) GO TO 4
      DO 9 I=1,N-1
         C = ABS(A(I,I))
         iRI = I
         iCI = I
         DO 10 J=I,N
            DO 11 K=I,N
               IF (ABS(A(J,K)) .GT. C) THEN
                  C = ABS(A(J,K))
                  iRI = J
                  iCI = K
               END IF
   11       CONTINUE
   10    CONTINUE
         IF (iRI .NE. I) THEN
            DET = -DET
            DO 13 J=1,N
               C = A(I,J)
               A(I,J) = A(iRI,J)
               A(iRI,J) = C
   13       CONTINUE
            IF (iCTR .GE. 0) THEN
               DO 14 J=1,M
                  C = B(I,J)
                  B(I,J) = B(iRI,J)
                  B(iRI,J) = C
   14          CONTINUE
            END IF
         END IF
         IF (iCI .NE. I) THEN
            DET = -DET
            DO 16 J=1,N
               C = A(J,I)
               A(J,I) = A(J,iCI)
               A(J,iCI) = C
   16       CONTINUE
         END IF
         iS(I) = iCI
         SUM = A(I,I)
         IF (ABS(SUM) .LE. EPS) then
            write(6,*) ' case 1. i,sum,eps', i,sum,eps
            GO TO 1
         endif
         DO 17 J=I+1,N
            C = A(J,I)/SUM
            DO 18 K=I+1,N
               A(J,K) = A(J,K)-C*A(I,K)
   18       CONTINUE
            IF (iCTR .GE. 0) THEN
               DO 19 K=1,M
                  B(J,K) = B(J,K)-C*B(I,K)
   19          CONTINUE
            END IF
   17    CONTINUE
    9 CONTINUE
    4 SUM = A(N,N)
      IF (ABS(SUM) .LE. EPS) then
            write(6,*) ' case 2. n,sum,eps', n,sum,eps
         GO TO 1
      endif
      IF (iCTR .GT. 0) GO TO 6
      DO 20 K=1,N
         DET = DET*A(K,K)
    2    IF (ABS(DET) .GT. 1.0D10) THEN
            DET = DET*1.0D-20
            iEX = iEX+20
            GO TO 2
         END IF
    3    IF (ABS(DET) .LE. 1.0D-10) THEN
            DET = DET*1.0D20
            iEX = iEX-20
            GO TO 3
         END IF
   20 CONTINUE
      IF (iCTR .LT. 0) GO TO 5
    6 DO 21 K=1,M
         B(N,K) = B(N,K)/SUM
   21 CONTINUE
      IF (N .EQ. 1) GO TO 5
      DO 22 I=N-1,1,-1
         C = A(I,I)
         DO 23 J=1,M
            SUM = B(I,J)
            DO 24 K=I+1,N
               SUM = SUM-A(I,K)*B(K,J)
   24       CONTINUE
            B(I,J) = SUM/C
   23    CONTINUE
   22 CONTINUE
      DO 25 K=N-1,1,-1
         I = iS(K)
         IF (I .NE. K) THEN
            DO 26 J=1,M
               C = B(I,J)
               B(I,J) = B(K,J)
               B(K,J) = C
   26       CONTINUE
         END IF
   25 CONTINUE
    5 iCTR = 0
      RETURN
    1 DET = 0.0D0
      iCTR = 1
      RETURN
    7 iCTR = -1
      RETURN
      END
