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

subroutine gen1overR3(Lhigh,oneoverR3)
!bs generates the radial integrals  for the one electron spin orbit integrals
!bs taken the 1/r**3 formula from the documentation and included additional
!bs factors for normalization

use AMFI_global, only: df, exponents, MxprimL, nprimit, Lmax
use Constants, only: Two, Quart, Pi
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: Lhigh
real(kind=wp), intent(out) :: oneoverR3(MxprimL*(MxprimL+1)/2,Lmax)
integer(kind=iwp) :: icount, iprim1, iprim2, L
real(kind=wp) :: alpha1, alpha2

do L=1,Lhigh
  icount = 0
  do iprim2=1,nprimit(L)
    alpha2 = exponents(iprim2,L)
    do iprim1=1,iprim2
      alpha1 = exponents(iprim1,L)
      icount = icount+1
      oneoverR3(icount,L) = sqrt(Two/Pi)*(df(L+L-2)*real(2**(L+3),kind=wp)*(alpha1*alpha2)**(Quart*real(L+L+3,kind=wp)))/ &
                            ((alpha1+alpha2)**L*df(L+L+1))
    end do
  end do
end do

return

end subroutine gen1overR3
