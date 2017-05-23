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
      CHARACTER*40 FUNCTION RDNAME (LUT,KEYW)
C
C...  reads a character item up wherever it is in a file in the form
C     ITEM_NAME  value
C
C     Returns ' ' when it fails to read the item up.
C
C       LUT       logical unit to be searched
C                 LUT>0 for formatted files, LUT<0 for unformatted files
C       KEYW         keyword
C
      INTEGER      LU,LUT
      CHARACTER*40 KEYW,VALUE,PIKNAM
      CHARACTER*80 LINE,BLANK
      LOGICAL      UNFORM
      EXTERNAL     PIKNAM
      BLANK =            '                                        '
      BLANK =BLANK(:40)//'                                        '
C
      IF (LUT.GT.0) THEN
        UNFORM=.FALSE.
        LU=LUT
      ELSE IF (LUT.LT.0) THEN
        UNFORM=.TRUE.
        LU=-LUT
C       if not IBM: no error messages for reading short records
      Else
        LU=0
        UNFORM=.TRUE.
        Write (6,*) 'RdName: LUT=0!'
        Call Abend
      ENDIF
C
      REWIND LU
1     CONTINUE
        IF (UNFORM) THEN
C...      unformatted files
          LINE=BLANK
          READ (LU,END=999)          LINE
          NOFB=MYLEN(LINE)
          IF (NOFB.EQ.0) NOFB=1
          IF (NOFB.LT.80) LINE=LINE(1:NOFB)
        ELSE
C...      formatted files
          READ (LU,'(A80)',END=999) LINE
        ENDIF
        VALUE=PIKNAM(LINE,KEYW)
        IF (VALUE.EQ.' ') GO TO 1
C
C...  the item has been found
      RDNAME=VALUE
      RETURN
C
C...  the item has not been found in the file
999   CONTINUE
      RDNAME=VALUE
      RETURN
C
      END
