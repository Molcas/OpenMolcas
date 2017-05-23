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
* Copyright (C) 1982,1983, Per Ake Malmqvist                           *
************************************************************************
      SUBROUTINE DOOL(NDIMEN,MDIM,N,M,A,B,DET,IPIV,JPIV,BUF)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NDIMEN,NDIMEN),B(NDIMEN,MDIM)
      DIMENSION IPIV(NDIMEN),JPIV(NDIMEN),BUF(NDIMEN)
C
C     SOLVES THE MATRIX EQUATION AX=B BY DOOLITTLE''S METHOD
C     ACTUAL DIMENSIONS ARE N*N AND N*M
C     ALLOCATED DIMENSIONS ARE NDIMEN*NDIMEN AND NDIMEN*MDIM
C     A AND B ARE DESTROYED, AND X IS RETURNED AS MATRIX B
C                                       (MALMQUIST 82-11-12)
C                                        (UPDATE 83-09-28)
C
C  EQUATION IS SOLVED BY FACTORIZING A=L*R IN SAME SPACE AS A.
C  PIVOTING IS ACHIEVED BY INDIRECT INDEXING.
C  FIRST PIVOTING INDICES ARE ASSIGNED START VALUES.
C
      IP=0 ! dummy initialize
      JP=0 ! dummy initialize
      DO I=1,N
       IPIV(I)=I
       JPIV(I)=I
      END DO
      DET=1.0D00
      DO I=1,N
C  NOW FIND BETTER PIVOT ELEMENT:
       AMAX=-1.0D00
       DO K=I,N
        DO L=I,N
        AM=ABS(A(IPIV(K),JPIV(L)))
        IF(AMAX.LE.AM) THEN
         AMAX=AM
         IP=K
         JP=L
        END IF
        END DO
       END DO
       IF(IP.NE.I) THEN
        DET=-DET
        IDUM=IPIV(I)
        IPIV(I)=IPIV(IP)
        IPIV(IP)=IDUM
       END IF
       IF(JP.NE.I) THEN
        DET=-DET
        IDUM=JPIV(I)
        JPIV(I)=JPIV(JP)
        JPIV(JP)=IDUM
       END IF
       IP=IPIV(I)
       JP=JPIV(I)
       DIAG=A(IP,JP)
       BUF(I)=DIAG
       DET=DET*DIAG
       DO K=I+1,N
        KP=IPIV(K)
        C=A(KP,JP)/DIAG
        A(KP,JP)=C
        DO L=I+1,N
         LP=JPIV(L)
         A(KP,LP)=A(KP,LP)-C*A(IP,LP)
        END DO
       END DO
      END DO
C
C  FIRST RESUBSTITUTION STEP:
C
      DO J=1,M
       DO I=2,N
        IP=IPIV(I)
        SUMMA=B(IP,J)
        DO K=1,I-1
         SUMMA=SUMMA-A(IP,JPIV(K))*B(IPIV(K),J)
        END DO
        B(IP,J)=SUMMA
       END DO
      END DO
C
C  SECOND RESUBSTITUTION STEP:
C
      DO J=1,M
       DO I=N,1,-1
        IP=IPIV(I)
        SUMMA=B(IP,J)
        DO K=I+1,N
         SUMMA=SUMMA-A(IP,JPIV(K))*B(IPIV(K),J)
        END DO
        B(IP,J)=SUMMA/BUF(I)
       END DO
      END DO
C
C  REORGANIZATION PART:
C
      DO J=1,M
       DO I=1,N
        BUF(I)=B(IPIV(I),J)
       END DO
       DO I=1,N
        B(JPIV(I),J)=BUF(I)
       END DO
      END DO
      RETURN
      END
