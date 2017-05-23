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
      SUBROUTINE CHO_OPEN(LUNIT,FNAME)
C
C     Purpose: open sequential unformatted fortran file.
C
      IMPLICIT NONE
      INTEGER  LUNIT
      CHARACTER*(*) FNAME

      INTEGER  LOCUNT, ISEED

      INTEGER ISFREEUNIT

      IF ((LUNIT.LT.1) .OR. (LUNIT.GT.99)) THEN
         LOCUNT = 7
      ELSE
         LOCUNT = LUNIT
      END IF

      ISEED  = LOCUNT
      LOCUNT = ISFREEUNIT(ISEED)
      call molcas_binaryopen_vanilla(Locunt, Fname)
c      OPEN(LOCUNT,FILE=FNAME,STATUS='UNKNOWN',FORM='UNFORMATTED')
      LUNIT = LOCUNT

      END
