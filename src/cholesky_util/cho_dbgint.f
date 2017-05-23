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
      SUBROUTINE CHO_DBGINT()
C
C     Purpose: regenerate and check integrals as specified in input.
C
      IMPLICIT NONE
#include "cholesky.fh"

      CALL QENTER('_DBGINT')

      IF (CHO_MINCHK) THEN ! minimal check
         CALL CHO_MCA_DBGINT_S(ICHKQ,NCHKQ,.TRUE.)
      ELSE ! check all (or the number of col. from input)
         CALL CHO_MCA_DBGINT_A()
      END IF

      CALL QEXIT('_DBGINT')

      END
