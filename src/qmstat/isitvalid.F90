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

subroutine IsItValid(Coo,CooRef,ValidOrNot)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
#include "maxi.fh"
real(kind=wp) :: Coo(3,MxCen), CooRef(MxCen,3)
logical(kind=iwp) :: ValidOrNot
integer(kind=iwp) :: i, j, k
real(kind=wp) :: dL_ref, dL_test
real(kind=wp), parameter :: dTroskel = 1.0e-4_wp

ValidOrNot = .true.
! Lengths.
outer: do i=1,4
  do j=i+1,5
    dL_test = Zero
    dL_ref = Zero
    do k=1,3
      dL_test = dL_test+(Coo(k,i)-Coo(k,j))**2
      dL_ref = dL_ref+(CooRef(i,k)-CooRef(j,k))**2
    end do
    if (abs(dL_test-dL_ref) > dTroskel) then
      ValidOrNot = .false.
      exit outer
    end if
  end do
end do outer

return

end subroutine IsItValid
