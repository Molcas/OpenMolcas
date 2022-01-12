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

subroutine Setup_OffAO()

use Basis_Info, only: dbsc, nCnttp, Shells
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iCnttp, kComp, kSh, lComp, lSh

! For some reason we need a counter which for a given shell index,
! kSh, counts how many angular functions have been in the
! proceeding shells. This is stored in kOffAO(kSh).
! Additionally, lOffAO gives the total number of angular functions
! in a given dbsc.

do iCnttp=1,nCnttp
  lComp = 0
  do lSh=0,dbsc(iCnttp)%nVal-1
    kSh = lSh+dbsc(iCnttp)%iVal
    if (Shells(kSh)%Prjct) then
      kComp = 2*lSh+1
    else
      kComp = (lSh+1)*(lSh+2)/2
    end if
    Shells(kSh)%kOffAO = lComp
    if ((Shells(kSh)%nBasis_C /= 0) .and. (Shells(kSh)%nExp /= 0)) lComp = lComp+kComp
  end do
  dbsc(iCnttp)%lOffAO = lComp
end do

return

end subroutine Setup_OffAO
