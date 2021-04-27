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

function qcomp(alpha,dfac,n,l,t,xk)
! Compute q(n,l) scaled by sqrt(pi)*exp(-t) to prevent overflows
! arguments are alpha, xk, and t=xk**2/(4*alpha)
! No restriction on the magnitude of t
! Increase dfac array to raise n, l restrictions

use Constants, only: Zero, One, Two, Four, Pi
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: qcomp
real(kind=wp), intent(in) :: alpha, dfac(*), t, xk
integer(kind=iwp), intent(in) :: n, l
integer(kind=iwp) :: num
real(kind=wp) :: dsum, fac1, fac2, fnum, prefac, term, xc, xden, xj, xkp, xnum
real(kind=wp), parameter :: accpow = 1.0e-14_wp, accasy = 1.0e-10_wp, a1rtpi = One/sqrt(Pi), &
                            tmin(0:8) = [31.0_wp,28.0_wp,25.0_wp,23.0_wp,22.0_wp,20.0_wp,19.0_wp,18.0_wp,15.0_wp]

if ((mod(n+l,2) == 0) .and. (n > l)) then

  ! use alternating series ((n+l <= 22) .and. (l <= 10))

  if (l == 0) then
    xkp = One
  else
    xkp = (xk/(alpha+alpha))**l
  end if
  prefac = xkp*dfac(n+l+1)/((alpha+alpha)**((n-l)/2)*sqrt(Four*alpha)*dfac(l+l+3))
  num = l-n+2
  xden = (l+l+3)
  term = One
  dsum = term
  xc = -One
  do
    if (num == 0) exit
    fnum = num
    term = term*fnum*t/(xden*xc)
    xc = xc-One
    dsum = dsum+term
    num = num+2
    xden = xden+Two
  end do

  qcomp = prefac*dsum

else if (t < tmin(min(n,8))) then

  ! use power series ((n+l <= 22) .and. (l <= 10))

  if (l == 0) then
    xkp = One
  else
    xkp = (xk/(alpha+alpha))**l
  end if
  prefac = exp(-t)*xkp/(alpha+alpha)**((n-l+1)/2)
  if (mod(n+l,2) == 0) then
    prefac = prefac/sqrt(Four*alpha)
  else
    prefac = a1rtpi*prefac
  end if
  xnum = (l+n-1)
  xden = (l+l+1)
  term = dfac(l+n+1)/dfac(l+l+3)
  dsum = term
  xj = Zero
  do
    xnum = xnum+Two
    xden = xden+Two
    xj = xj+One
    term = term*t*xnum/(xj*xden)
    dsum = dsum+term
    if ((term/dsum) <= accpow) exit
  end do
  qcomp = prefac*dsum

else

  ! use asymptotic series (arbitrary n, l)

  xkp = (xk/(alpha+alpha))**(n-2)
  prefac = xkp/((alpha+alpha)*sqrt(Four*alpha))
  dsum = One
  term = One
  fac1 = (l-n+2)
  fac2 = (1-l-n)
  xc = One
  do
    term = term*fac1*fac2/(Four*xc*t)
    if (term == Zero) exit
    dsum = dsum+term
    if (abs(term/dsum) < accasy) exit
    fac1 = fac1+Two
    fac2 = fac2+Two
    xc = xc+One
  end do
  qcomp = prefac*dsum

end if

return

end function qcomp
