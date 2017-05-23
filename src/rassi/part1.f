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
* Copyright (C) 1984,1989, Per Ake Malmqvist                           *
************************************************************************
C---------------------------------------------------------------
      SUBROUTINE PART1(NDIMEN,NBLOCK,NSIZE,SXY,B,A,SCR,IPIV,BUF)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SXY(NDIMEN,NDIMEN),A(NDIMEN,NDIMEN),B(NDIMEN,NDIMEN)
      DIMENSION SCR(NDIMEN,NDIMEN),BUF(NDIMEN),IPIV(NDIMEN,2)
      DIMENSION NSIZE(NBLOCK)
C
C  PURPOSE: SEE SUBROUTINE PART.
C  SUBDIVISION INTO TWO LEVELS OF ROUTINE CALLS IS MERELY TO
C  FACILITATE HANDLING OF SYMMETRY AND INDEXING.
C  ORIGINAL VERSION, MALMQUIST 84-04-04
C  RASSCF VERSION,   MALMQUIST 89-11-15
C---------------------------------------------------------------
C INITIALIZE A = INVERSE OF SXY, AND B = UNIT MATRIX:
      CALL QENTER('PART1')
      DO 20 I=1,NDIMEN
        DO 10 J=1,NDIMEN
          A(I,J)=0.0D00
          B(I,J)=0.0D00
          SCR(I,J)=SXY(I,J)
10      CONTINUE
        A(I,I)=1.0D00
        B(I,I)=1.0D00
20    CONTINUE
      CALL DOOL (NDIMEN,NDIMEN,NDIMEN,NDIMEN,SCR,A,DET,
     *           IPIV(1,1),IPIV(1,2),BUF)
C---------------------------------------------------------------
C LOOP BACKWARDS OVER THE BLOCKS. KEEP THREE LIMITS UPDATED:
C LIM1= POS IMMEDIATELY BEFORE CURRENT BLOCK, LIM2= BEGINNING
C OF CURRENT BLOCK, AND LIM3=END OF CURRENT BLOCK.
      LIM1=NDIMEN
      DO 100 K=NBLOCK,2,-1
        NSZ=NSIZE(K)
        LIM3=LIM1
        LIM1=LIM1-NSZ
        LIM2=LIM1+1
C---------------------------------------------------------------
C CALCULATE (INVERSE OF CURRENT A-BLOCK)*(ALL TO THE LEFT OF IT)
C AND PUT IT INTO B-MATRIX. THEN CLEAR ALL TO THE LEFT OF THE A-BLOCK.
        DO 50 I=LIM2,LIM3
          DO 30 J=LIM2,LIM3
            SCR(I,J)=A(I,J)
30        CONTINUE
          DO 40 J=1,LIM1
            B(I,J)=A(I,J)
            A(I,J)=0.0D00
40        CONTINUE
50      CONTINUE
        CALL DOOL(NDIMEN,NDIMEN,NSZ,LIM1,SCR(LIM2,LIM2),B(LIM2,1),DET,
     *            IPIV(1,1),IPIV(1,2),BUF)
C---------------------------------------------------------------
C THEN UPDATE THE COLUMNS OF A TO THE LEFT OF THE CURRENT BLOCK:
        DO 70 J=1,LIM1
          DO 70 I=1,LIM1
            T=A(I,J)
            DO 60 KK=LIM2,LIM3
              T=T-B(KK,J)*A(I,KK)
60          CONTINUE
            A(I,J)=T
70      CONTINUE
100   CONTINUE
C TRANSPOSE MATRIX B:
      DO 110 I=1,NDIMEN-1
        DO 110 J=I,NDIMEN
          T=B(I,J)
          B(I,J)=B(J,I)
          B(J,I)=T
110   CONTINUE
C---------------------------------------------------------------
C COMBINED LU-PARTITIONING AND UNITARY TRANSFORMATION OF A AND B:
      CALL LU2(NDIMEN,NBLOCK,NSIZE,B,A,BUF)
C
C     LU PARTITIONING OF THE MATRICES WAS DONE IN-PLACE.
C     NOW CHANGE SIGN OF LOWER-TRIANGULAR PARTS AND
C     INVERT UPPER-TRIANGULAR PARTS, AS INDICATED IN (VI.6):
C
      DO 200 I=2,NDIMEN
        DO 200 J=1,I-1
          A(I,J)=-A(I,J)
          B(I,J)=-B(I,J)
200   CONTINUE
      DO 240 L=NDIMEN,1,-1
        A(L,L)=1.0D0/A(L,L)
        B(L,L)=1.0D0/B(L,L)
        DO 210 M=L+1,NDIMEN
          A(L,M)=A(L,L)*A(L,M)
210       B(L,M)=B(L,L)*B(L,M)
        DO 230 K=1,L-1
          DO 220 M=L+1,NDIMEN
            A(K,M)=A(K,M)-A(K,L)*A(L,M)
220         B(K,M)=B(K,M)-B(K,L)*B(L,M)
          A(K,L)=-(A(L,L)*A(K,L))
          B(K,L)=-(B(L,L)*B(K,L))
230       CONTINUE
240     CONTINUE
      CALL QEXIT('PART1')
      RETURN
      END
