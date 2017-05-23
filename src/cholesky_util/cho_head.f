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
      SUBROUTINE CHO_HEAD(STRING,LINE,LENMAX,LUNIT)
C
C     Purpose: print a header.
C
      IMPLICIT NONE
      CHARACTER*(*) STRING
      CHARACTER*1   LINE
      INTEGER       LENMAX, LUNIT

      INTEGER LENSTR, LENTOT, I

      LENSTR = LEN(STRING)
      LENTOT = MIN(LENSTR,LENMAX-2)
      IF (LENTOT .GT. 0) THEN
         WRITE(LUNIT,'(//,2X,A)') STRING(1:LENTOT)
         WRITE(LUNIT,'(2X,80A)') (LINE,I=1,LENTOT)
      ELSE
         WRITE(LUNIT,'(//,2X,A,/)') STRING(1:LENSTR)
      END IF

      END
