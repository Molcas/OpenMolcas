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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine CASPT2_Grad_FwdCnt(iS,jS,kS,lS,LoadVec)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(inout) :: iS, jS, kS, lS
logical(kind=iwp), intent(out) :: LoadVec
integer(kind=iwp) :: iSold, kSold

LoadVec = .false.
iSold = iS
kSold = kS

lS = lS+1
if (iS == kS) then
  if (lS > jS) then
    jS = jS+1
    lS = 1
  end if
  if (jS > iS) then
    iS = iS+1
    jS = 1
    kS = 1
    lS = 1
  end if
else
  if (lS > kS) then
    jS = jS+1
    lS = 1
  end if
  if (jS > iS) then
    kS = kS+1
    jS = 1
    lS = 1
  end if
  if (kS > iS) then
    iS = iS+1
    jS = 1
    kS = 1
    lS = 1
  end if
end if

if ((iSold /= iS) .or. (kSold /= kS)) LoadVec = .true.

return

end subroutine CASPT2_Grad_FwdCnt
