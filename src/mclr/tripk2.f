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
      SUBROUTINE TRIPK2(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGN)
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
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AUTPAK(MATDIM,MATDIM),APAK(*)
C
      IF( IWAY .EQ. 1 ) THEN
        IJ = 0
        DO 100 J = 1,NDIM
          DO 50  I = J , NDIM
           APAK(IJ+I) = AUTPAK(I,J)
   50     CONTINUE
          IJ = IJ +NDIM-J
  100   CONTINUE
      END IF
C
      IF( IWAY .EQ. 2 ) THEN
        IJ = 0
        DO 200 J = 1,NDIM
          DO 150  I = J,NDIM
           AUTPAK(J,I) = SIGN*APAK(IJ+I)
           AUTPAK(I,J) = APAK(IJ+I)
  150     CONTINUE
          IJ = IJ + NDIM-J
  200   CONTINUE
      END IF
C
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' AUTPAK AND APAK FROM TRIPK2 '
        CALL WRTMAT(AUTPAK,NDIM,MATDIM,NDIM,MATDIM)
        CALL PRSM2(APAK,NDIM)
      END IF
C
      RETURN
      END
