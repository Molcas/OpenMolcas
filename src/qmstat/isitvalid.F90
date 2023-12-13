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

use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: Coo(3,5), CooRef(3,5)
logical(kind=iwp), intent(out) :: ValidOrNot
integer(kind=iwp) :: i, j
real(kind=wp) :: dL_ref, dL_test
real(kind=wp), parameter :: dTroskel = 1.0e-4_wp

ValidOrNot = .true.
! Lengths.
outer: do i=1,4
  do j=i+1,5
    dL_test = (Coo(1,i)-Coo(1,j))**2+(Coo(2,i)-Coo(2,j))**2+(Coo(3,i)-Coo(3,j))**2
    dL_ref = (CooRef(1,i)-CooRef(1,j))**2+(CooRef(2,i)-CooRef(2,j))**2+(CooRef(3,i)-CooRef(3,j))**2
    if (abs(dL_test-dL_ref) > dTroskel) then
      ValidOrNot = .false.
      exit outer
    end if
  end do
end do outer

return

end subroutine IsItValid
