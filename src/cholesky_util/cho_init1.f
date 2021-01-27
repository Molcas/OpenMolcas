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
      SUBROUTINE CHO_INIT1()
C
C     Purpose: initialize counter arrays.
C
      use ChoSwp, only: InfRed, InfVec
#include "implicit.fh"
#include "cholesky.fh"
#include "choglob.fh"
#include "cho_para_info.fh"

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CHO_INIT1')

      INTEGER  CHO_ISUMELM
      EXTERNAL CHO_ISUMELM

      IF (RSTCHO) THEN

C        Read restart information.
C        -------------------------

         CALL CHO_GETRSTC()
         NUMCHT = CHO_ISUMELM(NUMCHO,NSYM)

      ELSE

C        Initialize vector info and counters.
C        ------------------------------------

         CALL CHO_IZERO(INFVEC,SIZE(INFVEC))
         CALL CHO_IZERO(NUMCHO,NSYM)
         NUMCHT = 0

C        Initialize reduced set info.
C        ----------------------------

         CALL CHO_IZERO(INFRED,SIZE(INFRED))

C        Initialize global integral pass counter.
C        ----------------------------------------

         XNPASS = 0

      END IF

C     Parallel init.
C     --------------

      IF (Cho_Real_Par) CALL CHO_IZERO(MYNUMCHO,NSYM)

      END
