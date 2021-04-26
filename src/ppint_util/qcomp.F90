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

implicit real*8(a-h,o-z)
parameter(am1=-1.0d0,a0=0.0d0,accpow=1.0d-14,accasy=1.0d-10,a1rtpi=0.56418958354775629d0,a1=1.0d0,a2=2.0d0,a4=4.0d0)
dimension dfac(*)
dimension tmin(0:8)
data tmin/31.0d0,28.0d0,25.0d0,23.0d0,22.0d0,20.0d0,19.0d0,18.0d0,15.0d0/

if (mod(n+l,2) /= 0 .or. n <= l) go to 30

! use alternating series ((n+l <= 22) .and. (l <= 10))

if (l == 0) then
  xkp = a1
else
  xkp = (xk/(alpha+alpha))**l
end if
prefac = xkp*dfac(n+l+1)/((alpha+alpha)**((n-l)/2)*sqrt(a4*alpha)*dfac(l+l+3))
num = l-n+2
xden = (l+l+3)
term = a1
sum = term
xc = am1
10 continue
if (num /= 0) then
  fnum = num
  term = term*fnum*t/(xden*xc)
  xc = xc+am1
  sum = sum+term
  num = num+2
  xden = xden+a2
  go to 10
end if

qcomp = prefac*sum
return

30 continue
if (t < tmin(min(n,8))) go to 60

! use asymptotic series (arbitrary n, l)

xkp = (xk/(alpha+alpha))**(n-2)
prefac = xkp/((alpha+alpha)*sqrt(a4*alpha))
sum = a1
term = a1
fac1 = (l-n+2)
fac2 = (1-l-n)
xc = a1
40 continue
term = term*fac1*fac2/(a4*xc*t)
if (term == a0) go to 50
sum = sum+term
if (abs(term/sum) < accasy) go to 50
fac1 = fac1+a2
fac2 = fac2+a2
xc = xc+a1
go to 40
50 continue
qcomp = prefac*sum
return

! use power series ((n+l <= 22) .and. (l <= 10))

60 continue
if (l == 0) then
  xkp = a1
else
  xkp = (xk/(alpha+alpha))**l
end if
prefac = exp(-t)*xkp/(alpha+alpha)**((n-l+1)/2)
if (mod(n+l,2) == 0) then
  prefac = prefac/sqrt(a4*alpha)
else
  prefac = a1rtpi*prefac
end if
xnum = (l+n-1)
xden = (l+l+1)
term = dfac(l+n+1)/dfac(l+l+3)
sum = term
xj = a0
70 continue
xnum = xnum+a2
xden = xden+a2
xj = xj+a1
term = term*t*xnum/(xj*xden)
sum = sum+term
if ((term/sum) > accpow) go to 70
qcomp = prefac*sum

return

end function qcomp
