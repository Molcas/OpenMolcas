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

subroutine ZMatConv(LuWr,iAtom,iErr)

use ZMatConv_Mod, only: iZmat, Coords, Zmat
use Constants, only: Zero, One, deg2rad
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LuWR, iAtom
integer(kind=iwp), intent(out) :: iErr
integer(kind=iwp) :: i
real(kind=wp) :: arg, dCBiAtom, dCTiAtom, dSBiAtom, dSTiAtom, r, u1(3), u2(3), u3(3), u4(3), vj(3), vp(3)

iErr = 0

dCBiAtom = cos(Zmat(2,iAtom)*deg2rad)
dSBiAtom = sin(Zmat(2,iAtom)*deg2rad)
dCTiAtom = cos(Zmat(3,iAtom)*deg2rad)
dSTiAtom = sin(Zmat(3,iAtom)*deg2rad)
if (abs(dCBiAtom) < 1.0e-10_wp) dCBiAtom = Zero
if (abs(dSBiAtom) < 1.0e-10_wp) dSBiAtom = Zero
if (abs(dCTiAtom) < 1.0e-10_wp) dCTiAtom = Zero
if (abs(dSTiAtom) < 1.0e-10_wp) dSTiAtom = Zero

call vec(1.0e-6_wp,u1,iZmat(2,iAtom),iZmat(3,iAtom),iErr)
! Vet. u1(NB-NT)
if (iErr /= 0) then ! Vettore u1 nullo
  iErr = 1
  write(LuWr,*) ' [Z-Mat_Conv] Incipient floating point error detected for atom ',iAtom
  return
end if
call vec(1.0e-6_wp,u2,iZmat(1,iAtom),iZmat(2,iAtom),iErr)
! Vet. u2(NA,NB)
if (iErr /= 0) then ! Vettore u2 nullo
  iErr = 1
  write(LuWr,*) ' [Z-Mat_Conv] Incipient floating point error detected for atom ',iAtom
  return
end if
arg = One-(u1(1)*u2(1)+u1(2)*u2(2)+u1(3)*u2(3))**2
if (arg < Zero) then ! u1 e u2 allineati
  iErr = 1
  write(LuWr,*) ' [Z-Mat_Conv] Incipient floating point error detected for atom ',iAtom
  return
end if
r = sqrt(arg) ! (u1.u2)^0.5
if (r < 1.0e-6_wp) then ! r piccolo
  iErr = 1
  write(LuWr,*) ' [Z-Mat_Conv] Incipient floating point error detected for atom ',iAtom
  return
end if

call crprod(u1,u2,vp) ! Vettore piano perpendicolare u1 x u2
do i=1,3
  u3(i) = vp(i)/r
end do
call crprod(u3,u2,u4)
do i=1,3
  vj(i) = Zmat(1,iAtom)*(-u2(i)*dCBiAtom+u4(i)*dSBiAtom*dCTiAtom+u3(i)*dSBiAtom*dSTiAtom)
  Coords(i,iAtom) = Coords(i,iZmat(1,iAtom))+vj(i)
end do

return

end subroutine ZMatConv
