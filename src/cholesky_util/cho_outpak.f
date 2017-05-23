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
* Copyright (C) 1973, Nelson H. F. Beebe                               *
************************************************************************
      SUBROUTINE CHO_OUTPAK(AMATRX,NROW,NCTL,LUPRI)
C.......................................................................
C
C OUTPAK PRINTS A REAL SYMMETRIC MATRIX STORED IN ROW-PACKED LOWER
C
C TRIANGULAR FORM (SEE DIAGRAM BELOW) IN FORMATTED FORM WITH NUMBERED
C
C ROWS AND COLUMNS.  THE INPUT IS AS FOLLOWS:
C
C        AMATRX(*)...........PACKED MATRIX
C
C        NROW................NUMBER OF ROWS TO BE OUTPUT
C
C        NCTL................CARRIAGE CONTROL FLAG: 1 FOR SINGLE SPACE,
C                                                   2 FOR DOUBLE SPACE,
C                                                   3 FOR TRIPLE SPACE.
C
C THE MATRIX ELEMENTS ARE ARRANGED IN STORAGE AS FOLLOWS:
C
C        1
C        2    3
C        4    5    6
C        7    8    9   10
C       11   12   13   14   15
C       16   17   18   19   20   21
C       22   23   24   25   26   27   28
C       AND SO ON.
C
C OUTPAK IS SET UP TO HANDLE 6 COLUMNS/PAGE WITH A 6F20.14 FORMAT
C FOR THE COLUMNS.  IF A DIFFERENT NUMBER OF COLUMNS IS REQUIRED, CHANGE
C FORMATS 1000 AND 2000, AND INITIALIZE KCOL WITH THE NEW NUMBER OF
C COLUMNS.
C
C AUTHOR:  NELSON H.F. BEEBE, QUANTUM THEORY PROJECT, UNIVERSITY OF
C          FLORIDA, GAINESVILLE
C..........VERSION = 09/05/73/03
C.......................................................................
C
#include "implicit.fh"
      DIMENSION AMATRX(*)
      INTEGER BEGIN
      CHARACTER*1 ASA(3),BLANK,CTL
      CHARACTER   PFMT*20, COLUMN*8
      PARAMETER (ZERO=0.D00, KCOLP=4, KCOLN=6)
      PARAMETER (FFMIN=1.D-3, FFMAX = 1.D3)
      DATA COLUMN/'Column  '/, ASA/' ', '0', '-'/, BLANK/' '/
C
      IF (NCTL .LT. 0) THEN
         KCOL = KCOLN
      ELSE
         KCOL = KCOLP
      END IF
      MCTL = ABS(NCTL)
      IF ((MCTL.LE.3).AND.(MCTL.GT.0)) THEN
         CTL = ASA(MCTL)
      ELSE
         CTL = BLANK
      END IF
C
      J = NROW*(NROW+1)/2
      AMAX = ZERO
      DO 5 I=1,J
         AMAX = MAX( AMAX, ABS(AMATRX(I)) )
    5 CONTINUE
      IF (AMAX .EQ. ZERO) THEN
         WRITE (LUPRI,'(/T6,A)') 'Zero matrix.'
         GO TO 200
      END IF
      IF (FFMIN .LE. AMAX .AND. AMAX .LE. FFMAX) THEN
C        use F output format
         PFMT = '(A1,I7,2X,8F15.8)'
      ELSE
C        use 1PD output format
         PFMT = '(A1,I7,2X,1P,8D15.6)'
      END IF
C
C LAST IS THE LAST COLUMN NUMBER IN THE ROW CURRENTLY BEING PRINTED
C
      LAST = MIN(NROW,KCOL)
C
C BEGIN IS THE FIRST COLUMN NUMBER IN THE ROW CURRENTLY BEING PRINTED.
C
C.....BEGIN NON STANDARD DO LOOP.
      BEGIN= 1
 1050 NCOL = 1
         WRITE (LUPRI,1000) (COLUMN,I,I = BEGIN,LAST)
         DO 40 K = BEGIN,NROW
            KTOTAL = (K*(K-1))/2 + BEGIN - 1
            DO 10 I = 1,NCOL
               IF (AMATRX(KTOTAL+I) .NE. ZERO) GO TO 20
   10       CONTINUE
            GO TO 30
   20       WRITE (LUPRI,PFMT) CTL,K,(AMATRX(J+KTOTAL),J=1,NCOL)
   30       IF (K .LT. (BEGIN+KCOL-1)) NCOL = NCOL + 1
   40    CONTINUE
         LAST = MIN(LAST+KCOL,NROW)
         BEGIN= BEGIN + NCOL
      IF (BEGIN.LE.NROW) GO TO 1050
  200 CONTINUE
      RETURN
C
 1000 FORMAT (/12X,6(3X,A6,I4,2X),(3X,A6,I4))
C2000 FORMAT (A1,'Row',I4,2X,1P,8D15.6)
C2000 FORMAT (A1,I7,2X,1P,8D15.6)
      END
