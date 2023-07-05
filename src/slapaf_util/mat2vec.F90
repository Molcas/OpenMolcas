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

use Constants, only: Zero, One, Half, Pi
use Definitions, only: wp

implicit none
real(kind=wp) :: RotMat(3,3), RotVec(3), RotAng
real(kind=wp) :: C, cr, Pimth, sr, tn, tn2, trace, x1, x2, x3

trace = RotMat(1,1)+RotMat(2,2)+RotMat(3,3)
x1 = Half*(RotMat(3,2)-RotMat(2,3))
x2 = Half*(RotMat(1,3)-RotMat(3,1))
x3 = Half*(RotMat(2,1)-RotMat(1,2))
sr = sqrt(x1**2+x2**2+x3**2)
cr = Half*(trace-One)
if (0.05_wp*cr > sr) then
  tn = sr/cr
  tn2 = tn**2
  ! cos(theta) is positive, C is theta/sin(theta)
  C = (45045.0_wp-tn2*(15015.0_wp-tn2*(9009.0_wp-tn2*(6435.0_wp-tn2*(5005.0_wp-tn2*(4095.0_wp-tn2*3465.0_wp))))))/(45045.0_wp*cr)
  RotVec(1) = x1*C
  RotVec(2) = x2*C
  RotVec(3) = x3*C
  RotAng = sr*C
else if ((-0.05_wp*abs(cr) > sr) .and. (sr > Zero)) then
  tn = -sr/cr
  tn2 = tn**2
  ! cos(theta) is negative, sin(theta) is positive (always),
  ! Pimth is (Pi-theta)
  Pimth = tn*(45045.0_wp-tn2*(15015.0_wp-tn2*(9009.0_wp-tn2*(6435.0_wp-tn2*(5005.0_wp-tn2*(4095.0_wp-tn2*3465.0_wp))))))/45045.0_wp
  RotAng = Pi-Pimth
  RotVec(1) = (x1/sr)*RotAng
  RotVec(2) = (x2/sr)*RotAng
  RotVec(3) = (x3/sr)*RotAng
else if (sr /= Zero) then
  RotAng = atan2(sr,cr)
  RotVec(1) = RotAng*(x1/sr)
  RotVec(2) = RotAng*(x2/sr)
  RotVec(3) = RotAng*(x3/sr)
else
  RotVec(1) = Zero
  RotVec(2) = Zero
  RotVec(3) = Zero
  RotAng = Zero
end if

return

end subroutine mat2vec
