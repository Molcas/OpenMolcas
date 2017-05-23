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
* Copyright (C) 1984, Per Ake Malmqvist                                *
************************************************************************
C---------------------------------------------------------------
      SUBROUTINE LU2 (NDIMEN,NBLOCK,NSIZE,CXA,CYB,SCR)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CXA(NDIMEN,NDIMEN),CYB(NDIMEN,NDIMEN),SCR(NDIMEN)
      DIMENSION NSIZE(NBLOCK)
C  GIVES A SIMULTANEOUS LU-PARTITIONING OF MATRICES CXA,CYB IN THE
C  SENSE THAT CXA*X = L1*U1 AND CYB*X = L2*U2, WHERE X IS A BLOCK
C  UNITARY MATRIX. X IS NEVER FORMED, BUT AT EACH STEP SUCH A
C  TRANSFORMATION IS APPLIED THAT THE LU PARTITIONING CAN BE
C  APPLIED IN AN OPTIMAL WAY -- THIS IS REFERRED TO AS A
C  UNITARY PSEUDO-PIVOTATION. MATRICES CXA AND CYB ARE DESTROYED,
C  AND WILL CONTAIN THE NONTRIVIAL ELEMENTS OF THE TRIANGULAR
C  MATRICES.
C                                         ( MALMQUIST 84-01-16 )
      NTOT=0
      DO IBLOCK=1,NBLOCK
        NTOT=NTOT+NSIZE(IBLOCK)
      END DO
      IF(NTOT.NE.NDIMEN) WRITE(6,*)' ERROR: NTOT.NE.NDIMEN IN LU2.'
      IEND=0
      DO IBLOCK=1,NBLOCK
        ISTA=IEND+1
        IEND=IEND+NSIZE(IBLOCK)
        DO II=ISTA,IEND
          S1=0.0D00
          S2=0.0D00
          S3=0.0D00
          DO J=II,IEND
            S1=S1+CXA(II,J)**2
            S2=S2+CYB(II,J)**2
            S3=S3+CXA(II,J)*CYB(II,J)
          END DO
          THR=1.0D-06
          IF((S1.LT.THR).OR.(S2.LT.THR)) THEN
            WRITE(6,*)' RASSI CANNOT CONTINUE. THE PROBLEM AT HAND'
            WRITE(6,*)' IS PROBABLY NOT SOLUBLE. THE TWO ORBITAL'
            WRITE(6,*)' SPACES ARE TOO DISSIMILAR.'
            WRITE(6,*)' LU PARTITIONING IS TROUBLESOME. DIAGONAL'
            WRITE(6,*)' ELEMENT NR.',II,' IS TOO SMALL:'
            WRITE(6,*)' IN MATRIX CXA IT IS', SQRT(S1)
            WRITE(6,*)' IN MATRIX CYB IT IS', SQRT(S2)
            WRITE(6,*)' EVEN AFTER OPTIMAL PIVOT-TRANSFORMATION.'
            CALL ABEND()
          END IF
          X1=1.0D00/SQRT(S1)
          X2=1.0D00/SIGN(SQRT(S2),S3)
          DO I=II,IEND
            SCR(I)=X1*CXA(II,I)+X2*CYB(II,I)
          END DO
          S=2.0D0*(1.0D00+X1*X2*S3)
          X=1.0D00/SIGN(SQRT(S),SCR(II))
          DO I=II,IEND
            SCR(I)=X*SCR(I)
          END DO
          X=1.0D00/(1.0D00+SCR(II))
          DO I=1,IEND
            S=0.0D00
            DO J=II,IEND
              S=S+CXA(I,J)*SCR(J)
            END DO
            S2=0.0D00
            DO J=II+1,IEND
              S2=S2+CXA(I,J)*SCR(J)
            END DO
            S2=CXA(I,II)+X*S2
            CXA(I,II)=S
            DO J=II+1,IEND
              CXA(I,J)=CXA(I,J)-S2*SCR(J)
            END DO
          END DO
          DO I=1,IEND
            S=0.0D00
            DO J=II,IEND
              S=S+CYB(I,J)*SCR(J)
            END DO
            S2=0.0D00
            DO J=II+1,IEND
              S2=S2+CYB(I,J)*SCR(J)
            END DO
            S2=CYB(I,II)+X*S2
            CYB(I,II)=S
            DO J=II+1,IEND
            CYB(I,J)=CYB(I,J)-S2*SCR(J)
            END DO
          END DO
          X=1.0D00/CXA(II,II)
          DO I=II+1,IEND
            CXA(I,II)=X*CXA(I,II)
            DO J=II+1,NTOT
              CXA(I,J)=CXA(I,J)-CXA(I,II)*CXA(II,J)
            END DO
          END DO
          X=1.0D00/CYB(II,II)
          DO I=II+1,IEND
            CYB(I,II)=X*CYB(I,II)
            DO J=II+1,NTOT
              CYB(I,J)=CYB(I,J)-CYB(I,II)*CYB(II,J)
            END DO
          END DO
        END DO
      END DO
      RETURN
      END
