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
* Copyright (C) 1994, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE GATRMT(MATIN,NROWIN,NCOLIN,MATUT,NROWUT,NCOLUT,
     &                  ISCA,SSCA)
*
* Gather rows of transposed matrix MATIN  to  MATUT
*
* MATUT(I,J) = SSCA(I)*MATIN(J,ISCA(I)),(ISCA(I) .ne. 0 )
*
*
* Jeppe Olsen, Getting LUCIA in shape , Feb1994
*
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 MATIN,MATUT
*.Input
      DIMENSION ISCA(*),SSCA(*),MATIN(NCOLIN,NROWIN)
*. ( MATIN Transposed )
*.Output
      DIMENSION MATUT(NROWUT,NCOLUT)
*
*. To get rid of annoying and incorrect compiler warnings
      ICOFF = 0
C     LBLK = 100
      LBLK = 40
      NBLK = NCOLUT/LBLK
      IF(LBLK*NBLK.LT.NCOLUT) NBLK = NBLK + 1
      DO ICBL = 1, NBLK
        IF(ICBL.EQ.1) THEN
          ICOFF = 1
        ELSE
          ICOFF = ICOFF + LBLK
        END IF
        ICEND = MIN(ICOFF+LBLK-1,NCOLUT)
        DO I = 1, NROWUT
          IF(ISCA(I).NE.0) THEN
            S = SSCA(I)
            IROW = ISCA(I)
            DO J = ICOFF,ICEND
              MATUT(I,J) = S*MATIN(J,IROW)
            END DO
          ELSE IF (ISCA(I).EQ.0) THEN
            DO J = ICOFF,ICEND
              MATUT(I,J) = 0.0D0
            END DO
          END IF
        END DO
      END DO
*
      RETURN
      END
      SUBROUTINE SCARMT(MATIN,NROWIN,NCOLIN,MATUT,NROWUT,NCOLUT,
     &                  ISCA,SSCA)
*
* Scatter-add  rows of MATIN to transposed matrix MATUT

*  MATUT(J,ISCA(I)) = MATUT(J,ISCA(I)) + SSCA(I)*MATIN(I,J)
*  ( if INDEX(I).ne.0 )
*
* Jeppe Olsen, Getting LUCIA in shape , Feb1994
*
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 MATIN,MATUT
*.Input
      DIMENSION ISCA(*),SSCA(*),MATIN(NROWIN,NCOLIN)
*.Input and Output
      DIMENSION MATUT(NCOLUT,NROWUT)
*.                 (MATUT transposed !)
*. To get rid of annoying and incorrect compiler warnings
      ICINOF = 0
*
C     LBLK = 100
      LBLK = 40
      NBLK = NCOLIN/LBLK
      IF(LBLK*NBLK.LT.NCOLIN) NBLK = NBLK + 1
      DO ICINBL = 1, NBLK
        IF(ICINBL.EQ.1) THEN
          ICINOF = 1
        ELSE
          ICINOF = ICINOF + LBLK
        END IF
        ICINEN = MIN(ICINOF+LBLK-1,NCOLIN)
        DO I = 1, NROWIN
          IF(ISCA(I).NE.0) THEN
            S = SSCA(I)
            IROW = ISCA(I)
            DO ICOL = ICINOF,ICINEN
              MATUT(ICOL,IROW) = MATUT(ICOL,IROW)+S*MATIN(I,ICOL)
            END DO
          END IF
        END DO
      END DO
*
      RETURN
      END
