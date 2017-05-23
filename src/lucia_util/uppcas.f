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
      SUBROUTINE UPPCAS(LINE,LENGTH)
*
* Convert letters in character string LINE to upper case
*
* very stupid and not vectorized !
*
      CHARACTER*(*) LINE
      PARAMETER (NCHAR = 41)
      CHARACTER*1 LOWER(NCHAR)
      CHARACTER*1 UPPER(NCHAR)
*
      DATA LOWER/'a','b','c','d','e',
     &           'f','g','h','i','j',
     &           'k','l','m','n','o',
     &           'p','q','r','s','t',
     &           'u','v','w','x','y',
     &           'z','+','-','<','>',
     &           '=','0','1','2','3',
     &           '4','5','6','7','8',
     &           '9'/
      DATA UPPER/'A','B','C','D','E',
     &           'F','G','H','I','J',
     &           'K','L','M','N','O',
     &           'P','Q','R','S','T',
     &           'U','V','W','X','Y',
     &           'Z','+','-','<','>',
     &           '=','0','1','2','3',
     &           '4','5','6','7','8',
     &           '9'/
*
      DO 100 ICHA = 1, LENGTH
        DO 50 I = 1,NCHAR
          IF(LINE(ICHA:ICHA).EQ.LOWER(I))
     &    LINE(ICHA:ICHA) = UPPER(I)
   50   CONTINUE
  100 CONTINUE
*
      RETURN
      END
