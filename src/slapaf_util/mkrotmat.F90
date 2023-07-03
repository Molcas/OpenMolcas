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

implicit real*8(a-h,o-z)
real*8 RotMat(3,3)
real*8 RotVec(3)

!call RecPrt('mkRotMat: RotVec',' ',RotVec,1,3)

Q = RotVec(1)**2+RotVec(2)**2+RotVec(3)**2
XN = sqrt(Q)
if (Q < 0.01d0) then
  c0 = 1.d0-(Q/2.d0)*(1.d0-(Q/12.d0)*(1.d0-(Q/30.d0)*(1.d0-(Q/56.d0))))
  c1 = 1.d0-(Q/6.d0)*(1.d0-(Q/20.d0)*(1.d0-(Q/42.d0)*(1.d0-(Q/72.d0))))
  c2 = (1.d0-(Q/12.d0)*(1.d0-(Q/30.d0)*(1.d0-(Q/56.d0)*(1.d0-(Q/90.d0)))))/2.d0
else
  cs = cos(XN)
  sn = sin(XN)
  c0 = cs
  c1 = sn/XN
  c2 = (1.0d0-cs)/(XN**2)
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
do i=1,3
  do j=1,3
    RotMat(i,j) = RotMat(i,j)+c2*RotVec(i)*RotVec(j)
  end do
end do
do i=1,3
  do j=1,3
    sum = 0.0d0
    if (i == j) sum = -1.0d0
    do k=1,3
      sum = sum+RotMat(i,k)*RotMat(j,k)
    end do
    if (abs(sum) > 1.0D-10) then
      call WarningMessage(2,'Error in RotDer')
      write(6,*) ' MKROTMAT: ON check sum error=',sum
      call Abend()
    end if
  end do
end do

return

end subroutine mkRotMat
