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

function Index_NoSym(iCntr,iCmp,iCnt,iAng,iR,cIndex,iBas,nBas)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: Index_NoSym
integer(kind=iwp), intent(in) :: iCntr, iCmp, iCnt, iAng, iR, nBas
integer(kind=iwp), intent(inout) :: cIndex(5,nBas), iBas
integer(kind=iwp) :: i

outer: do
  Index_NoSym = 0
  do i=1,iBas
    if ((cIndex(1,i) == iCntr) .and. (cIndex(2,i) == iCmp) .and. (cIndex(3,i) == iCnt) .and. (cIndex(4,i) == iAng) .and. &
        (cIndex(5,i) == iR)) then
      Index_NoSym = i
      exit outer
    end if
  end do

  iBas = iBas+1
  cIndex(1,iBas) = iCntr
  cIndex(2,iBas) = iCmp
  cIndex(3,iBas) = iCnt
  cIndex(4,iBas) = iAng
  cIndex(5,iBas) = iR
end do outer

end function Index_NoSym
