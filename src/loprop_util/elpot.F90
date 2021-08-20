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

real*8 function ElPot(r,rinv,x,y,z,dMullig,lMax,A,chP,lDOrNot1,lDOrNot2)
! The electric potential with diffuse s- and p-functions. No
! higher than d-functions (with non-zero trace).

implicit real*8(a-h,o-z)

#include "WrkSpc.fh"
#include "warnings.h"

dimension A(2), dMullig((lMax*(lMax**2+6*lMax+11)+6)/6)
dimension dL2(6), dL3(10), dL4(15), dL5(21)

logical lDOrNot1, lDOrNot2

ElPot = 0.0d0
if (lMax >= 0) then
  if (lDOrNot1) then
    dCh = chP*rinv
    dCh = dCh+dMullig(1)*rinv*(1.0d0-(1.0d0+A(1)*r)*exp(-2.0d0*A(1)*r))
    ElPot = dCh
  else
    dL2(1) = dMullig(1)+chP
    ElPot = ElPot+ElPointPot(rinv,x,y,z,0,dL2)
  end if
end if
if (lMax >= 1) then
  if (lDOrNot2) then
    ar = A(2)*r
    dM = (x*dMullig(2)+y*dMullig(3)+z*dMullig(4))*rinv**3*(1.0d0-(1.0d0+2.0d0*ar+2.0d0*ar**2+ar**3)*exp(-2.0d0*ar))
    ElPot = ElPot+dM
  else
    dL2(1) = dMullig(2)
    dL2(2) = dMullig(3)
    dL2(3) = dMullig(4)
    ElPot = ElPot+ElPointPot(rinv,x,y,z,1,dL2)
  end if
end if
if (lMax >= 2) then
  dL2(1) = dMullig(5)
  dL2(2) = dMullig(6)
  dL2(3) = dMullig(7)
  dL2(4) = dMullig(8)
  dL2(5) = dMullig(9)
  dL2(6) = dMullig(10)
  ElPot = ElPot+ElPointPot(rinv,x,y,z,2,dL2)
end if
if (lMax >= 3) then
  dL3(1) = dMullig(11)
  dL3(2) = dMullig(12)
  dL3(3) = dMullig(13)
  dL3(4) = dMullig(14)
  dL3(5) = dMullig(15)
  dL3(6) = dMullig(16)
  dL3(7) = dMullig(17)
  dL3(8) = dMullig(18)
  dL3(9) = dMullig(19)
  dL3(10) = dMullig(20)
  ElPot = ElPot+ElPointPot(rinv,x,y,z,3,dL3)
end if
if (lMax >= 4) then
  dL4(1) = dMullig(21)
  dL4(2) = dMullig(22)
  dL4(3) = dMullig(23)
  dL4(4) = dMullig(24)
  dL4(5) = dMullig(25)
  dL4(6) = dMullig(26)
  dL4(7) = dMullig(27)
  dL4(8) = dMullig(28)
  dL4(9) = dMullig(29)
  dL4(10) = dMullig(30)
  dL4(11) = dMullig(31)
  dL4(12) = dMullig(32)
  dL4(13) = dMullig(33)
  dL4(14) = dMullig(34)
  dL4(15) = dMullig(35)
  ElPot = ElPot+ElpointPot(rinv,x,y,z,4,dL4)
end if
if (lMax >= 5) then
  dL5(1) = dMullig(36)
  dL5(2) = dMullig(37)
  dL5(3) = dMullig(38)
  dL5(4) = dMullig(39)
  dL5(5) = dMullig(40)
  dL5(6) = dMullig(41)
  dL5(7) = dMullig(42)
  dL5(8) = dMullig(43)
  dL5(9) = dMullig(44)
  dL5(10) = dMullig(45)
  dL5(11) = dMullig(46)
  dL5(12) = dMullig(47)
  dL5(13) = dMullig(48)
  dL5(14) = dMullig(49)
  dL5(15) = dMullig(50)
  dL5(16) = dMullig(51)
  dL5(17) = dMullig(52)
  dL5(18) = dMullig(53)
  dL5(19) = dMullig(54)
  dL5(20) = dMullig(55)
  dL5(21) = dMullig(56)
  ElPot = ElPot+ElpointPot(rinv,x,y,z,5,dL5)
end if
if (lMax >= 6) then
  write(6,*)
  write(6,*) 'Oops! You hit the roof with respect to angular momentum. Lower that, or do some programming.'
  call Quit(_RC_GENERAL_ERROR_)
end if

return

end function ElPot
