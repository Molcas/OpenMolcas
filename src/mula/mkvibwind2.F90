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
! Copyright (C) 2009, Giovanni Ghigo                                   *
!***********************************************************************

subroutine MkVibWind2(nYes,iMaxYes,max_nOrd,lVec,VibWind2)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nYes, max_nOrd, lVec(0:max_nOrd)
integer(kind=iwp), intent(out) :: iMaxYes, VibWind2(nYes)
integer(kind=iwp) :: iOrd, iYes

iMaxYes = 0
iYes = 1
do iOrd=0,max_nOrd
  if (lVec(iOrd) == 1) then
    VibWind2(iYes) = iOrd
    iYes = iYes+1
    iMaxYes = max(iMaxYes,iOrd)
  end if
end do

return

end subroutine MkVibWind2
