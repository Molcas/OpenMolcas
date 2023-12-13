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

subroutine CHO_RSTD_GETIND2()
!
! Purpose: read mapping arrays for diagonal restart.

use Cholesky, only: IndRed, IndRSh, LuRed, nnBstRT, nnShl, nSym
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: IADR, IOPT, LEN0, LEN1, LEN2

LEN0 = NSYM*NNSHL
LEN1 = NNBSTRT(1)
LEN2 = LEN1
IOPT = 2
IADR = LEN0
call IDAFILE(LURED,IOPT,INDRED,LEN1,IADR)
IOPT = 2
IADR = LEN0+LEN1
call IDAFILE(LURED,IOPT,INDRSH,LEN2,IADR)

end subroutine CHO_RSTD_GETIND2
