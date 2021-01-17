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
      SUBROUTINE CHO_RSTD_GETIND1()
C
C     Purpose: read and set some index arrays for diagonal restart.
C
      use ChoSwp, only: nnBstRSh
#include "implicit.fh"
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

C     Read info.
C     ----------

      IOPT = 2
      IADR = 0
      LTOT = NSYM*NNSHL
      CALL IDAFILE(LURED,IOPT,NNBSTRSH,LTOT,IADR)

C     Set up IIBSTRSH etc.
C     --------------------

      IRED = 1
      CALL CHO_SETREDIND(IWORK(ip_IIBSTRSH),NNBSTRSH,
     &                   NSYM,NNSHL,IRED)

      END
