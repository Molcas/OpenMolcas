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

subroutine CHO_RESETCNF()
!
! Purpose: reset configuration of decomposition to that read from
!          restart file. Original configuration will be saved in
!          restart common block.

use Cholesky, only: DAMP, SCDIAG, SPAN, THRCOM, THRDIAG, THRNEG, TOONEG, WARNEG, XDAMP, XSCDIAG, XSPAN, XTHRCOM, XTHRDIAG, &
                    XTHRNEG, XTOONEG, XWARNEG
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: Tmp(8)
logical(kind=iwp) :: lTmp

Tmp(1) = XTHRCOM
Tmp(2) = XTHRDIAG
Tmp(3:4) = XDAMP(:)
Tmp(5) = XSPAN
Tmp(6) = XTHRNEG
Tmp(7) = XWARNEG
Tmp(8) = XTOONEG
lTmp = XSCDIAG

XTHRCOM = THRCOM
XTHRDIAG = THRDIAG
XDAMP(:) = DAMP(:)
XSPAN = SPAN
XTHRNEG = THRNEG
XWARNEG = WARNEG
XTOONEG = TOONEG
XSCDIAG = SCDIAG

THRCOM = Tmp(1)
THRDIAG = Tmp(2)
DAMP(:) = Tmp(3:4)
SPAN = Tmp(5)
THRNEG = Tmp(6)
WARNEG = Tmp(7)
TOONEG = Tmp(8)
SCDIAG = lTmp

end subroutine CHO_RESETCNF
