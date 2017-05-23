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
      SUBROUTINE TRIPK3(  AUTPAK,    APAK,    IWAY,  MATDIM,    NDIM,
     &                      SIGN)
C
C
C.. REFORMATING BETWEEN LOWER TRIANGULAR PACKING
C   AND FULL MATRIX FORM FOR A SYMMETRIC OR ANTI SYMMETRIC MATRIX
C
C   IWAY = 1 : FULL TO PACKED
C              LOWER HALF OF AUTPAK IS STORED IN APAK
C   IWAY = 2 : PACKED TO FULL FORM
C              APAK STORED IN LOWER HALF
C               SIGN * APAK TRANSPOSED IS STORED IN UPPPER PART
C.. NOTE : COLUMN WISE STORAGE SCHEME IS USED FOR PACKED BLOCKS
*
* Some considerations on cache minimization used for IMET = 2 Loop
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AUTPAK(MATDIM,MATDIM),APAK(*)
*. To get rid of annoying and incorrect compiler warnings
      IOFF = 0
      JOFF = 0
*
*. Packing : No problem with cache misses
*
      IF( IWAY .EQ. 1 ) THEN
        IJ = 0
        DO J = 1,NDIM
          DO I = J , NDIM
            APAK(IJ+I) = AUTPAK(I,J)
          END DO
          IJ = IJ +NDIM-J
        END DO
      END IF
*
* Unpacking : cache misses can occur so two routes
*
      IF( IWAY .EQ. 2 ) THEN
*. Use blocked algorithm
      IMET = 2
      IF(IMET.EQ.1) THEN
*. No blocking
        IJ = 0
        DO J = 1,NDIM
          DO I = J,NDIM
           AUTPAK(J,I) = SIGN*APAK(IJ+I)
           AUTPAK(I,J) = APAK(IJ+I)
          END DO
          IJ = IJ + NDIM-J
        END DO
      ELSE IF (IMET .EQ. 2 ) THEN
*. Blocking
        LBLK = 40
        NBLK = MATDIM/LBLK
        IF(LBLK*NBLK.LT.MATDIM) NBLK = NBLK + 1
        DO JBLK = 1, NBLK
          IF(JBLK.EQ.1) THEN
            JOFF = 1
          ELSE
            JOFF = JOFF + LBLK
          END IF
          JEND = MIN(JOFF+LBLK-1,MATDIM)
          DO IBLK = JBLK, NBLK
            IF(IBLK.EQ.JBLK) THEN
              IOFF = JOFF
            ELSE
              IOFF = IOFF + LBLK
            END IF
            IEND = MIN(IOFF+LBLK-1,MATDIM)
              DO J = JOFF,JEND
                IF(IBLK.EQ.JBLK) THEN
                  IOFF2 = J
                ELSE
                  IOFF2 = IOFF
                END IF
                IJOFF = (J-1)*MATDIM-J*(J-1)/2
                DO I = IOFF2,IEND
                  AUTPAK(J,I) = SIGN*APAK(IJOFF+I)
                  AUTPAK(I,J) = APAK(IJOFF+I)
                END DO
              END DO
*. End of loop over I and J
            END DO
          END DO
*. End of loop over blocks of I and J
        END IF
      END IF
*
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' AUTPAK AND APAK FROM TRIPK3 '
        CALL WRTMAT(AUTPAK,NDIM,MATDIM,NDIM,MATDIM)
        CALL PRSM2(APAK,NDIM)
      END IF
*
      RETURN
      END
