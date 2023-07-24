!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      SUBROUTINE CHO_RSTD_GETIND1()
!
!     Purpose: read and set some index arrays for diagonal restart.
!
      use ChoSwp, only: nnBstRSh
      Implicit Real*8 (a-h,o-z)
#include "cholesky.fh"

!     Read info.
!     ----------

      IOPT = 2
      IADR = 0
      LTOT = NSYM*NNSHL
      CALL IDAFILE(LURED,IOPT,NNBSTRSH,LTOT,IADR)

!     Set up IIBSTRSH etc.
!     --------------------

      IRED = 1
      CALL CHO_SETREDIND(IRED)

      END
