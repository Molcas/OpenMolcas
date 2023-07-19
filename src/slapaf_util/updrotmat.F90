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
real(kind=wp), intent(in) :: SmallRot(3)
real(kind=wp), intent(inout) :: RotMat(3,3)
integer(kind=iwp) :: i, j, k
real(kind=wp) :: rsum, SmallMat(3,3), tmp(3,3)

call mkRotMat(SmallRot,SmallMat)
tmp(:,:) = Zero
do i=1,3
  do j=1,3
    do k=1,3
      tmp(i,j) = tmp(i,j)+RotMat(i,k)*SmallMat(k,j)
    end do
  end do
end do
RotMat(:,:) = tmp(:,:)
! Check for orthonormality:
do i=1,3
  do j=1,3
    if (i == j) then
      rsum = -One
    else
      rsum = Zero
    end if
    do k=1,3
      rsum = rsum+RotMat(i,k)*RotMat(j,k)
    end do
    if (abs(rsum) > 1.0e-10_wp) write(u6,*) ' UPDROTMAT ON check sum error:',rsum
  end do
end do

return

end subroutine updRotMat
