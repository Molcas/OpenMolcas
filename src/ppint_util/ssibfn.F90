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

subroutine ssibfn(nmax,x,ssi)
! scaled spherical i Bessel functions

use Constants, only: Zero, One, Two, Three, OneHalf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nmax
real(kind=wp), intent(in) :: x
real(kind=wp), intent(out) :: ssi(nmax+1)
integer(kind=iwp) :: n
real(kind=wp) :: ak, akpkm2, akqkm2, aprod, ex, f2kp1, f2kp3, f2nm1, f2np1, f2np3, pk, pkm1, qk, qkm1, x2, xmin

x2 = x*x
xmin = (abs(3*nmax-1))

if (x <= xmin) then

  n = nmax
  f2np1 = (n+n+1)
  f2kp3 = f2np1
  pkm1 = Zero
  pk = One
  qkm1 = One
  qk = One
  aprod = One
  do
    f2kp1 = f2kp3
    f2kp3 = f2kp3+Two
    ak = x2/(f2kp1*f2kp3)
    akpkm2 = ak*pkm1
    pkm1 = pk
    pk = pkm1+akpkm2
    akqkm2 = ak*qkm1
    qkm1 = qk
    qk = qkm1+akqkm2
    aprod = ak*aprod
    if (((pk*qkm1)+aprod) == (pk*qkm1)) exit
  end do
  ssi(n+1) = pk/qk
  do
    if (n == 0) exit
    n = n-1
    f2np3 = f2np1
    f2np1 = f2np1-Two
    ssi(n+1) = (f2np1*f2np3)/((f2np1*f2np3)+x2*ssi(n+2))
  end do
  ssi(1) = ssi(1)/(One+x*ssi(1))
  do n=1,nmax
    ssi(n+1) = ssi(n+1)*ssi(n)
  end do

else

  if (x >= 20.0_wp) then
    ex = Zero
  else
    ex = exp(-(x+x))
  end if
  ssi(1) = (One-ex)/(x+x)
  if (nmax == 0) return
  ssi(2) = OneHalf*(One+ex+(ex-One)/x)/x2
  if (nmax == 1) return
  f2np1 = Three
  do n=2,nmax
    f2nm1 = f2np1
    f2np1 = f2np1+Two
    ssi(n+1) = (ssi(n-1)-ssi(n))*(f2nm1*f2np1)/x2
  end do

end if

return

end subroutine ssibfn
