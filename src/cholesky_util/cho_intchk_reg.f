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
      SUBROUTINE CHO_INTCHK_REG(LABEL,ISHLCD,ISHLAB)
C
C     Purpose: register a shell quadruple (CD|AB) for minimal integral
C              check using LABEL to keep track of its origin.
C
      use ChoArr, only: iSP2F
#include "implicit.fh"
      CHARACTER*8 LABEL
#include "cholesky.fh"

      CHARACTER*14 SECNAM
      PARAMETER (SECNAM = 'CHO_INTCHK_REG')

C     Check shell pair index.
C     -----------------------

      IF (ISHLCD.LT.1 .OR. ISHLCD.GT.NNSHL) THEN
         CALL CHO_QUIT('Shell index error 1 in '//SECNAM,103)
      END IF
      IF (ISHLAB.LT.1 .OR. ISHLAB.GT.NNSHL) THEN
         CALL CHO_QUIT('Shell index error 2 in '//SECNAM,103)
      END IF

C     Registration.
C     -------------

      CALL CHO_INVPCK(ISP2F(ISHLCD),ISHLC,ISHLD,.TRUE.)
      CALL CHO_INVPCK(ISP2F(ISHLAB),ISHLA,ISHLB,.TRUE.)
      CALL CHO_INTCHK_ID_OF(LABEL,ID,1)
      IF (ID.LT.1 .OR. ID.GT.NCHKQ) THEN
         ID = NCHKQ + 1 ! junk yard
      END IF
      ICHKQ(1,ID) = ISHLC
      ICHKQ(2,ID) = ISHLD
      ICHKQ(3,ID) = ISHLA
      ICHKQ(4,ID) = ISHLB

      END
