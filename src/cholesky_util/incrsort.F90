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
! Copyright (C) 2016, Giovanni Li Manni                                *
!***********************************************************************

subroutine IncrSort(EVal,EVec,dimens)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimens
real(kind=wp), intent(inout) :: EVal(dimens), EVec(dimens,dimens)
integer(kind=iwp) :: i, k, j, l
real(kind=wp) :: swap

do i=1,dimens-1
  k = i
  do j=i+1,dimens
    if (EVal(j) > EVal(k)) k = j
  end do
  if (k /= i) then
    Swap = EVal(k)
    EVal(k) = EVal(i)
    EVal(i) = Swap
    do l=1,dimens
      Swap = EVec(l,k)
      EVec(L,K) = EVec(l,i)
      EVec(L,I) = Swap
    end do
  end if
end do

return

end subroutine IncrSort
