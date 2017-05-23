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
      SUBROUTINE CHO_CLOSE(LUNIT,STAT)
C
C     Purpose: close sequential unformatted fortran file.
C
      IMPLICIT NONE
      INTEGER  LUNIT
      CHARACTER*(*) STAT

      IF ((LUNIT.LT.1) .OR. (LUNIT.GT.99)) THEN
         CALL CHO_QUIT('CHO_CLOSE: unit out of bounds!',104)
      ELSE
         CLOSE(LUNIT,STATUS=STAT)
         LUNIT = -1
      END IF

      END
