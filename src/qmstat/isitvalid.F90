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

implicit real*8(a-h,o-z)
#include "maxi.fh"
parameter(dTroskel=1d-4)
dimension Coo(MxCen,3), CooRef(MxCen,3)
logical ValidOrNot

ValidOrNot = .true.
! Lengths.
outer: do i=1,4
  do j=i+1,5
    dL_test = 0.0d0
    dL_ref = 0.0d0
    do k=1,3
      dL_test = dL_test+(Coo(i,k)-Coo(j,k))**2
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
