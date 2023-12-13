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

subroutine mkRotMat(rotvec,rotmat)

use Constants, only: Zero, One, Two, Six, Twelve, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: RotVec(3)
real(kind=wp), intent(out) :: RotMat(3,3)
integer(kind=iwp) :: i, j, k
real(kind=wp) :: c0, c1, c2, cs, Q, rsum, sn, XN

!call RecPrt('mkRotMat: RotVec',' ',RotVec,1,3)

Q = RotVec(1)**2+RotVec(2)**2+RotVec(3)**2
XN = sqrt(Q)
if (Q < 0.01_wp) then
  c0 = One-(Q/Two)*(One-(Q/Twelve)*(One-(Q/30.0_wp)*(One-(Q/56.0_wp))))
  c1 = One-(Q/Six)*(One-(Q/20.0_wp)*(One-(Q/42.0_wp)*(One-(Q/72.0_wp))))
  c2 = (One-(Q/Twelve)*(One-(Q/30.0_wp)*(One-(Q/56.0_wp)*(One-(Q/90.0_wp)))))*Half
else
  cs = cos(XN)
  sn = sin(XN)
  c0 = cs
  c1 = sn/XN
  c2 = (One-cs)/(XN**2)
end if
RotMat(1,1) = c0
RotMat(2,2) = c0
RotMat(3,3) = c0
RotMat(3,2) = c1*RotVec(1)
RotMat(1,3) = c1*RotVec(2)
RotMat(2,1) = c1*RotVec(3)
RotMat(2,3) = -c1*RotVec(1)
RotMat(3,1) = -c1*RotVec(2)
RotMat(1,2) = -c1*RotVec(3)
do j=1,3
  RotMat(:,j) = RotMat(:,j)+c2*RotVec(:)*RotVec(j)
end do
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
    if (abs(rsum) > 1.0e-10_wp) then
      call WarningMessage(2,'Error in mkRotMat')
      write(6,*) ' MKROTMAT: ON check sum error=',rsum
      call Abend()
    end if
  end do
end do

return

end subroutine mkRotMat
