!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2022, Jie J. Bao                                       *
!***********************************************************************
!*****************************************************************
! history:                                                       *
! Jie J. Bao, on Apr. 12, 2022, created this file.               *
!*****************************************************************

subroutine CalcQaa(Qaa,DDg,lRoots,nDDg)

use Constants, only: Zero, Half

implicit none
integer lRoots, nDDg
real*8 DDg(nDDg)
real*8 Qaa
integer iState, iLoc, Int1, lRoots2

lRoots2 = lRoots**2
Int1 = (lRoots2+1)*(lRoots+1)
Qaa = Zero
do iState=1,lRoots
  iLoc = (iState-1)*Int1+1
  Qaa = Qaa+DDg(iLoc)
end do
Qaa = Qaa*Half

return

end subroutine
