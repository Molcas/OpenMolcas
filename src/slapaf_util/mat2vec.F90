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

subroutine mat2vec(RotMat,RotVec,RotAng)

implicit real*8(a-h,o-z)
parameter(Pi=3.14159265358979323846d0)
dimension RotMat(3,3), RotVec(3)

trace = RotMat(1,1)+RotMat(2,2)+RotMat(3,3)
x1 = 0.5d0*(RotMat(3,2)-RotMat(2,3))
x2 = 0.5d0*(RotMat(1,3)-RotMat(3,1))
x3 = 0.5d0*(RotMat(2,1)-RotMat(1,2))
sr = sqrt(x1**2+x2**2+x3**2)
cr = 0.5d0*(trace-1.0d0)
if (0.05d0*cr > sr) then
  tn = sr/cr
  tn2 = tn**2
  ! cos(theta) is positive, C is theta/sin(theta)
  C = (45045.0d0-tn2*(15015.0d0-tn2*(9009.0d0-tn2*(6435.0d0-tn2*(5005.0d0-tn2*(4095.0d0-tn2*3465.0d0))))))/(45045.0d0*cr)
  RotVec(1) = x1*C
  RotVec(2) = x2*C
  RotVec(3) = x3*C
  RotAng = sr*C
else if ((-0.05d0*abs(cr) > sr) .and. (sr > 0.0d0)) then
  tn = -sr/cr
  tn2 = tn**2
  ! cos(theta) is negative, sin(theta) is positive (always),
  ! Pimth is (Pi-theta)
  Pimth = tn*(45045.0d0-tn2*(15015.0d0-tn2*(9009.0d0-tn2*(6435.0d0-tn2*(5005.0d0-tn2*(4095.0d0-tn2*3465.0d0))))))/45045.0d0
  RotAng = Pi-Pimth
  RotVec(1) = (x1/sr)*RotAng
  RotVec(2) = (x2/sr)*RotAng
  RotVec(3) = (x3/sr)*RotAng
else if (sr /= 0.0d0) then
  RotAng = atan2(sr,cr)
  RotVec(1) = RotAng*(x1/sr)
  RotVec(2) = RotAng*(x2/sr)
  RotVec(3) = RotAng*(x3/sr)
else
  RotVec(1) = 0.0d0
  RotVec(2) = 0.0d0
  RotVec(3) = 0.0d0
  RotAng = 0.0d0
end if

return

end subroutine mat2vec
