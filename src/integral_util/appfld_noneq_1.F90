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
subroutine AppFld_NonEq_1(Cavxyz,radius,Eps,lmax,EpsInf)

use Index_Functions, only: nTri3_Elem1
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lMax
real(kind=wp), intent(inout) :: Cavxyz(nTri3_Elem1(lMax))
real(kind=wp), intent(in) :: Radius, Eps, EpsInf
integer(kind=iwp) :: ip, l
real(kind=wp) :: Fact, rInv, rPoti
real(kind=wp), allocatable :: CavSph(:)
real(kind=wp), external :: DblFac

#ifdef _DEBUGPRINT_
call RecPrt('Multipole Moments',' ',Cavxyz,nTri3_Elem1(lMax),1)
#endif

! Backtransform from cartesian to spherical harmonics

call mma_allocate(CavSph,(lmax+1)**2,Label='CavSph')
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
  fact = f(Eps,l)*(One-f(EpsInf,l)/f(Eps,l))**2
  rpoti = rinv*fact*DblFac(2*l-1)
  CavSph(ip:ip+2*l) = rpoti*CavSph(ip:ip+2*l)
  ip = ip+2*l+1
end do

! Transform electric field components from spherical harmonics
! to cartesians.

call Tranca(Cavxyz,Cavsph,lmax,.false.)
call mma_deallocate(Cavsph)

#ifdef _DEBUGPRINT_
call RecPrt('Electric Field',' ',Cavxyz,nTri3_Elem1(lMax),1)
#endif

return

contains

pure function f(Eps_,l)

  real(kind=wp) :: f
  real(kind=wp), intent(in) :: Eps_
  integer(kind=iwp), intent(in) :: l

  f = (real(1+l,kind=wp)*(Eps_-One))/(real(1+l,kind=wp)*Eps_+real(l,kind=wp))

end function f

end subroutine AppFld_NonEq_1
