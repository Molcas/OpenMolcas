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
subroutine AppFld(Cavxyz,radius,Eps,lmax,EpsInf,NonEq)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lMax
real(kind=wp), intent(inout) :: Cavxyz((lMax+1)*(lMax+2)*(lMax+3)/6)
real(kind=wp), intent(in) :: Radius, Eps, EpsInf
logical(kind=iwp), intent(in) :: NonEq
integer(kind=iwp) :: ip, l
real(kind=wp) :: Fact, RInv, rPoti
real(kind=wp), allocatable :: CavSph(:)
real(kind=wp), external :: DblFac
! Statement function
real(kind=wp) :: F, Eps_
f(Eps_,l) = (real(1+l,kind=wp)*(Eps_-One))/(real(1+l,kind=wp)*Eps_+real(l,kind=wp))

#ifdef _DEBUGPRINT_
call RecPrt('Multipole Moments',' ',Cavxyz,(lMax+1)*(lMax+2)*(lMax+3)/6,1)
#endif

! Backtransform from cartesian to spherical harmonics

call mma_allocate(CavSph,(lmax+1)**2,Label='CavSph')
call Tranca(Cavxyz,Cavsph,lmax,.true.)
#ifdef _DEBUGPRINT_
call RecPrt(' CavSph',' ',Cavsph,(lMax+1),1)
#endif

! Evaluate the electric field components at the origin.
! This is identical to the charge distribution on the
! boundary of the cavity!

ip = 1
if (NonEq) then
  do l=0,lmax
    rinv = One/radius**(2*l+1)
    fact = (Two*f(EpsInf,l)-f(EpsInf,l)**2/f(Eps,l))
    rpoti = rinv*fact*DblFac(2*l-1)
    call DScal_(2*l+1,rpoti,Cavsph(ip),1)
    ip = ip+2*l+1
  end do
else
  do l=0,lmax
    rinv = One/radius**(2*l+1)
    fact = f(Eps,l)
    rpoti = rinv*fact*DblFac(2*l-1)
    call DScal_(2*l+1,rpoti,Cavsph(ip),1)
    ip = ip+2*l+1
  end do
end if

! Transform electric field components from spherical harmonics
! to cartesians.

call Tranca(Cavxyz,Cavsph,lmax,.false.)
call mma_deallocate(CavSph)

#ifdef _DEBUGPRINT_
call RecPrt('Electric Field',' ',Cavxyz,(lMax+1)*(lMax+2)*(lMax+3)/6,1)
#endif

return

end subroutine AppFld
