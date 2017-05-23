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
      SUBROUTINE CHO_MCA_CALCINT(ISHLAB)
C
C     Purpose: calculate qualified integral columns from
C              shell pair distribution (**|ISHLA ISHLB).
C
      IMPLICIT NONE
      INTEGER ISHLAB
#include "cholesky.fh"

      CHARACTER*15 SECNAM
      PARAMETER (SECNAM = 'CHO_MCA_CALCINT')

      IF (IFCSEW .EQ. 1) THEN ! store full shell quadruple
         CALL CHO_MCA_CALCINT_1(ISHLAB)
      ELSE IF (IFCSEW .EQ. 2) THEN ! store only reduced sets
         CALL CHO_MCA_CALCINT_2(ISHLAB)
      ELSE ! this is an error
         CALL CHO_QUIT('IFCSEW out of bounds in '//SECNAM,105)
      END IF

      END
