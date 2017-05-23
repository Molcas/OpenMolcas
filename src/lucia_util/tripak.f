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
      SUBROUTINE TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM)
C
C ( NOT A SIMPLIFIED VERSION OF TETRAPAK )
C
C.. REFORMATING BETWEEN LOWER TRIANGULAR PACKING
C   AND FULL MATRIX FORM FOR A SYMMETRIC MATRIX
C
C   IWAY = 1 : FULL TO PACKED
C   IWAY = 2 : PACKED TO FULL FORM
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AUTPAK(MATDIM,MATDIM),APAK(*)
C
      IF( IWAY .EQ. 1 ) THEN
        IJ = 0
        DO 100 I = 1,NDIM
          DO 50  J = 1, I
           APAK(IJ+J) = AUTPAK(J,I)
   50     CONTINUE
          IJ = IJ + I
  100   CONTINUE
      END IF
C
      IF( IWAY .EQ. 2 ) THEN
        IJ = 0
        DO 200 I = 1,NDIM
          DO 150  J = 1, I
           AUTPAK(I,J) = APAK(IJ+J)
           AUTPAK(J,I) = APAK(IJ+J)
  150     CONTINUE
          IJ = IJ + I
  200   CONTINUE
      END IF
C
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' AUTPAK AND APAK FROM TRIPAK '
        CALL WRTMAT(AUTPAK,NDIM,MATDIM,NDIM,MATDIM)
        CALL PRSYM(APAK,NDIM)
      END IF
C
      RETURN
      END
