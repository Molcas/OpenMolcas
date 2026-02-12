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
      use definitions, only: iwp, wp, u6
      use constants, only: Zero, One, Two
      IMPLICIT NONE
      integer(kind=iwp), intent(in):: NDIMEN,NBLOCK
      integer(kind=iwp), intent(in):: NSIZE(NBLOCK)
      real(kind=wp), intent(inout):: CXA(NDIMEN,NDIMEN),
     &                               CYB(NDIMEN,NDIMEN)
      real(kind=wp), intent(inout):: SCR(NDIMEN)

      integer(kind=iwp) NTOT,IBLOCK,IEND,ISTA,II,J,I
      real(kind=wp) S1,S2,S3,X1,X2,S,X
      real(kind=wp), parameter:: THR=1.0D-06
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
          S1=Zero
          S2=Zero
          S3=Zero
          DO J=II,IEND
            S1=S1+CXA(II,J)**2
            S2=S2+CYB(II,J)**2
            S3=S3+CXA(II,J)*CYB(II,J)
          END DO
          IF((S1.LT.THR).OR.(S2.LT.THR)) THEN
            WRITE(u6,*)' RASSI CANNOT CONTINUE. THE PROBLEM AT HAND'
            WRITE(u6,*)' IS PROBABLY NOT SOLUBLE. THE TWO ORBITAL'
            WRITE(u6,*)' SPACES ARE TOO DISSIMILAR.'
            WRITE(u6,*)' LU PARTITIONING IS TROUBLESOME. DIAGONAL'
            WRITE(u6,*)' ELEMENT NR.',II,' IS TOO SMALL:'
            WRITE(u6,*)' IN MATRIX CXA IT IS', SQRT(S1)
            WRITE(u6,*)' IN MATRIX CYB IT IS', SQRT(S2)
            WRITE(u6,*)' EVEN AFTER OPTIMAL PIVOT-TRANSFORMATION.'
            CALL ABEND()
          END IF
          X1=One/SQRT(S1)
          X2=One/SIGN(SQRT(S2),S3)
          DO I=II,IEND
            SCR(I)=X1*CXA(II,I)+X2*CYB(II,I)
          END DO
          S=Two*(One+X1*X2*S3)
          X=One/SIGN(SQRT(S),SCR(II))
          DO I=II,IEND
            SCR(I)=X*SCR(I)
          END DO
          X=One/(One+SCR(II))
          DO I=1,IEND
            S=Zero
            DO J=II,IEND
              S=S+CXA(I,J)*SCR(J)
            END DO
            S2=Zero
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
            S=Zero
            DO J=II,IEND
              S=S+CYB(I,J)*SCR(J)
            END DO
            S2=Zero
            DO J=II+1,IEND
              S2=S2+CYB(I,J)*SCR(J)
            END DO
            S2=CYB(I,II)+X*S2
            CYB(I,II)=S
            DO J=II+1,IEND
            CYB(I,J)=CYB(I,J)-S2*SCR(J)
            END DO
          END DO
          X=One/CXA(II,II)
          DO I=II+1,IEND
            CXA(I,II)=X*CXA(I,II)
            DO J=II+1,NTOT
              CXA(I,J)=CXA(I,J)-CXA(I,II)*CXA(II,J)
            END DO
          END DO
          X=One/CYB(II,II)
          DO I=II+1,IEND
            CYB(I,II)=X*CYB(I,II)
            DO J=II+1,NTOT
              CYB(I,J)=CYB(I,J)-CYB(I,II)*CYB(II,J)
            END DO
          END DO
        END DO
      END DO

      END SUBROUTINE LU2
