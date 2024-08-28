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

!#define _DEBUGPRINT_
subroutine AppFld_2(Cavxyz,Cavsph,radius,Eps,lmax,EpsInf,NonEq)

use Constants, only: One
use Definitions, only: wp

implicit none
integer lMax
real*8 Cavxyz((lMax+1)*(lMax+2)*(lMax+3)/6), Cavsph((lMax+1)**2)
real*8 Radius, Eps, EpsInf
logical NonEq
integer ip, l
real*8 rInv, Fact, rPoti, DblFac, F
! Statement function
f(Eps,l) = (real(1+l,kind=wp)*(Eps-One))/(real(1+l,wp)*Eps+real(l,wp))

#ifdef _DEBUGPRINT_
call RecPrt('Multipole Moments',' ',Cavxyz,(lMax+1)*(lMax+2)*(lMax+3)/6,1)
#endif

! Backtransform from cartesian to spherical harmonics

call Tranca(Cavxyz,Cavsph,lmax,.true.)
#ifdef _DEBUGPRINT_
call RecPrt(' CavSph',' ',Cavsph,(lMax+1)**2,1)
#endif

! Evaluate the electric field components at the origin.
! This is identical to the charge distribution on the
! boundary of the cavity!

ip = 1
do l=0,lmax
  rinv = One/radius**(2*l+1)
  fact = F(Eps,l)-F(EpsInf,l)-(F(EpsInf,l)-F(EpsInf,l)**2/F(Eps,l))
  rpoti = rinv*fact*DblFac(2*l-1)
  call DScal_(2*l+1,rpoti,Cavsph(ip),1)
  ip = ip+2*l+1
end do

! Transform electric field components from spherical harmonics
! to cartesians.

call Tranca(Cavxyz,Cavsph,lmax,.false.)

#ifdef _DEBUGPRINT_
call RecPrt('Electric Field',' ',Cavxyz,(lMax+1)*(lMax+2)*(lMax+3)/6,1)
#endif

return
! Avoid unused argument warnings
if (.false.) call Unused_logical(NonEq)

end subroutine AppFld_2
