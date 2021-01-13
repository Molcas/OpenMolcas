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
      CHARACTER*40 FUNCTION PIKNAM (LINED,KEYW)
C
C...  picks a character item up from a line like
C     ITEM_NAME  value  COMMENTS
C
C     returns ' ' when it fails to pick the item up.
C
C       LINED        input line to be searched
C       KEYW         keyword
C
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*40 KEYW
      CHARACTER*80 LINE,LINED
      LOGICAL FOUND
C
      IK=0
      LL=0
      IL=0
C
      LINE=LINED
      LK=1
      DO 10 I=40,1,-1
        IF (KEYW(I:I).NE.' ') THEN
          LK=I
          GO TO 11
        ENDIF
10    CONTINUE
11    CONTINUE
      DO 12 I=1,40
        IF (KEYW(I:I).NE.' ') THEN
          IK=I
          GO TO 13
        ENDIF
12    CONTINUE
13    CONTINUE
C
      ISTRT=INDEX(LINE,KEYW(IK:LK))
      IF (ISTRT.EQ.0) THEN
C...    the keyword has not been found in the line.
        PIKNAM=' '
        RETURN
      ENDIF
      ISTRT=ISTRT+LK-IK+1
      IF (LINE(ISTRT:ISTRT).NE.' ') THEN
C...    the keyword is not followed by a blank character.
        PIKNAM=' '
        RETURN
      ELSE IF (ISTRT.GE.80) THEN
C...    the line finishes after the keyword.
        PIKNAM=' '
        RETURN
      ELSE
        ISTRT=ISTRT+1
      ENDIF
C
      FOUND=.FALSE.
      DO 20 I=ISTRT,80
        IF (FOUND) THEN
          IF (LINE(I:I).EQ.' ') THEN
            LL=I-1
            GO TO 21
          ENDIF
        ELSE
          IF (LINE(I:I).NE.' ') THEN
            FOUND=.TRUE.
            IL=I
          ENDIF
        ENDIF
20    CONTINUE
21    CONTINUE
C
      IF (FOUND) THEN
        PIKNAM=LINE(IL:LL)
          WRITE (LINE,601) LK,LL-IL+1
C         WRITE (6,LINE) KEYW(:LK),PIKNAM
      ELSE
        PIKNAM=' '
      ENDIF
C
      RETURN
601   FORMAT ('('' PIKNAM:         '',A',I2.2,',1X,A',I2.2,')')
      END
