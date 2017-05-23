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
      SUBROUTINE CHO_RSTD_GETIND2()
C
C     Purpose: read mapping arrays for diagonal restart.
C
#include "implicit.fh"
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      LEN0 = NSYM*NNSHL
      LEN1 = NNBSTRT(1)
      LEN2 = LEN1
      IOPT = 2
      IADR = LEN0
      CALL IDAFILE(LURED,IOPT,IWORK(ip_INDRED),LEN1,IADR)
      IOPT = 2
      IADR = LEN0 + LEN1
      CALL IDAFILE(LURED,IOPT,IWORK(ip_INDRSH),LEN2,IADR)

      END
