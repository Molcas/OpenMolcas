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
! Copyright (C) 2007,2008, Roland Lindh                                *
!***********************************************************************

subroutine Remove_High_Exponents(iD,nD,List2,mData,nTheta_All)
!***********************************************************************
!                                                                      *
!     Experimental code to be used with care.                          *
!                                                                      *
!***********************************************************************

use Basis_Info, only: Shells
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(inout) :: nD, iD(nD)
integer(kind=iwp), intent(in) :: mData, nTheta_All, List2(mData,nTheta_All)
integer(kind=iwp) :: i, iTheta_All, k, kAng, kShll, l, lAng, lShll, mD
logical(kind=iwp) :: Skip

call iVcPrt('Remove_High_Exponents: iD',' ',iD,nD)
mD = nD
i = 1
do
  iTheta_All = iD(i)
  Skip = .false.
  kAng = List2(1,iTheta_All)
  lAng = List2(2,iTheta_All)
  k = List2(5,iTheta_All)
  l = List2(6,iTheta_All)
  kShll = List2(7,iTheta_All)
  lShll = List2(8,iTheta_All)
  if (kAng == lAng) then
    l = List2(6,iTheta_All)
    Skip = ((k == 1) .and. (l == 1)) .and. (Shells(kShll)%nExp /= 1)
  else
    Skip = (l == 1) .and. (Shells(lShll)%nExp /= 1)
  end if
  if (Skip) then
    if (mD == i) then
      mD = mD-1
      exit
    end if
    iD(i:mD-1) = iD(i+1:mD)
    mD = mD-1
    cycle
  end if
  i = i+1
  if (i > mD) exit
end do
nD = mD
call iVcPrt('Remove_High_Exponents: iD',' ',iD,nD)

return

end subroutine Remove_High_Exponents
