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

implicit real*8(a-h,o-z)
implicit integer(i-n)
dimension u1(3), u2(3), u3(3), u4(3), vj(3), vp(3)
#include "g_zmatconv.fh"

iErr = 0

torad = 3.1415926536d0/180.0d0
dCBiAtom = cos(torad*Zmat(iAtom,2))
dSBiAtom = sin(torad*Zmat(iAtom,2))
dCTiAtom = cos(torad*Zmat(iAtom,3))
dSTiAtom = sin(torad*Zmat(iAtom,3))
if (abs(dCBiAtom) < 1.0d-10) dCBiAtom = 0.0d0
if (abs(dSBiAtom) < 1.0d-10) dSBiAtom = 0.0d0
if (abs(dCTiAtom) < 1.0d-10) dCTiAtom = 0.0d0
if (abs(dSTiAtom) < 1.0d-10) dSTiAtom = 0.0d0

call vec(1.0d-06,u1,iZmat(iAtom,2),iZmat(iAtom,3),iErr)
! Vet. u1(NB-NT)
if (iErr /= 0) goto 9990 ! Vettore u1 nullo
call vec(1.0d-06,u2,iZmat(iAtom,1),iZmat(iAtom,2),iErr)
! Vet. u2(NA,NB)
if (iErr /= 0) goto 9990 ! Vettore u2 nullo
arg = 1.0d0-(u1(1)*u2(1)+u1(2)*u2(2)+u1(3)*u2(3))**2
if (arg < 0.0d0) goto 9990 ! u1 e u2 allineati
r = sqrt(arg) ! (u1.u2)^0.5
if (r < 1.0d-6) goto 9990 ! r piccolo

call crprod(u1,u2,vp) ! Vettore piano perpendicolare u1 x u2
do i=1,3
  u3(i) = vp(i)/r
end do
call crprod(u3,u2,u4)
do i=1,3
  vj(i) = Zmat(iAtom,1)*(-u2(i)*dCBiAtom+u4(i)*dSBiAtom*dCTiAtom+u3(i)*dSBiAtom*dSTiAtom)
  Coords(iAtom,i) = Coords(iZmat(iAtom,1),i)+vj(i)
end do

goto 9999

9990 iErr = 1
write(LuWr,*) ' [Z-Mat_Conv] Incipient floating point error detected for atom ',iAtom

9999 return

end subroutine ZMatConv
