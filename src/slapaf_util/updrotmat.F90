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

subroutine updRotMat(SmallRot,RotMat)

use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: SmallRot(3), RotMat(3,3)
integer(kind=iwp) :: i, j, k
real(kind=wp) :: rsum, SmallMat(3,3), tmp(3,3)

call mkRotMat(SmallRot,SmallMat)
do i=1,3
  do j=1,3
    rsum = Zero
    do k=1,3
      rsum = rsum+RotMat(i,k)*SmallMat(k,j)
    end do
    tmp(i,j) = rsum
  end do
end do
do i=1,3
  do j=1,3
    RotMat(i,j) = tmp(i,j)
  end do
end do
! Check for orthonormality:
do i=1,3
  do j=1,3
    rsum = Zero
    if (i == j) rsum = -One
    do k=1,3
      rsum = rsum+RotMat(i,k)*RotMat(j,k)
    end do
    if (abs(rsum) > 1.0e-10_wp) then
      write(u6,*) ' UPDROTMAT ON check sum error:',rsum
    end if
  end do
end do

return

end subroutine updRotMat
