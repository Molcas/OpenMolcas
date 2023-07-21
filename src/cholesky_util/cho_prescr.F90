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

subroutine Cho_PreScr(Thr1,Thr2)
!
! Purpose: read integral prescreening thresholds from common block.

use Gateway_Info, only: ThrInt, CutInt

real*8 Thr1, Thr2

Thr1 = CutInt
Thr2 = ThrInt

end subroutine Cho_PreScr
