!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      CHARACTER*40 FUNCTION RDNAME (LUT,KEYW)
!
!...  reads a character item up wherever it is in a file in the form
!     ITEM_NAME  value
!
!     Returns ' ' when it fails to read the item up.
!
!       LUT       logical unit to be searched
!                 LUT>0 for formatted files, LUT<0 for unformatted files
!       KEYW         keyword
!
      INTEGER      LU,LUT
      CHARACTER*40 KEYW,VALUE,PIKNAM
      CHARACTER*80 LINE,BLANK
      LOGICAL      UNFORM
      EXTERNAL     PIKNAM
      BLANK =            '                                        '
      BLANK =BLANK(:40)//'                                        '
!
      IF (LUT.GT.0) THEN
        UNFORM=.FALSE.
        LU=LUT
      ELSE IF (LUT.LT.0) THEN
        UNFORM=.TRUE.
        LU=-LUT
!       if not IBM: no error messages for reading short records
      Else
        LU=0
        UNFORM=.TRUE.
        Write (6,*) 'RdName: LUT=0!'
        Call Abend
      ENDIF
!
      REWIND LU
1     CONTINUE
        IF (UNFORM) THEN
!...      unformatted files
          LINE=BLANK
          READ (LU,END=999)          LINE
        ELSE
!...      formatted files
          READ (LU,'(A80)',END=999) LINE
        ENDIF
        VALUE=PIKNAM(LINE,KEYW)
        IF (VALUE.EQ.' ') GO TO 1
!
!...  the item has been found
      RDNAME=VALUE
      RETURN
!
!...  the item has not been found in the file
999   CONTINUE
      RDNAME=VALUE
      RETURN
!
      END
