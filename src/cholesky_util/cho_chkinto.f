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
      SUBROUTINE CHO_CHKINTO(XINT,DIAG,ISYM,NERR,TOL,REPORT)
C
C     Purpose: check diagonals in qualified integral columns against
C              original diagonal (read in here).
C
#include "implicit.fh"
      DIMENSION XINT(*), DIAG(*)
      LOGICAL   REPORT
#include "cholesky.fh"

      IOPT = 2
      CALL CHO_IODIAG(DIAG,IOPT)
      CALL CHO_P_CHKINT(XINT,DIAG,ISYM,NERR,TOL,REPORT)

      END
