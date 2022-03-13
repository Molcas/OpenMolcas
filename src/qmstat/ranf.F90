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

function ranf(idum)

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: ranf
integer(kind=iwp) :: idum
integer(kind=iwp), parameter :: IA1 = 40014, IA2 = 40692, IM1 = 2147483563, IM2 = 2147483399, IMM1 = IM1-1, IQ1 = 53668, &
                                IQ2 = 52774, IR1 = 12211, IR2 = 3791, NTAB = 32, NDIV = 1 +int(real(IMM1,kind=wp)/NTAB)
integer(kind=iwp) :: idum2 = 123456789, iv(NTAB) = 0, iy = 0, j, k
real(kind=wp), parameter :: AM = One/real(IM1,kind=wp), EPS = 1.2e-7_wp, RNMX = One-EPS

if (idum <= 0) then
  idum = max(-idum,1)
  idum2 = idum
  do j=NTAB+8,1,-1
    k = idum/IQ1
    idum = IA1*(idum-k*IQ1)-k*IR1
    if (idum < 0) idum = idum+IM1
    if (j <= NTAB) iv(j) = idum
  end do
  iy = iv(1)
end if
k = idum/IQ1
idum = IA1*(idum-k*IQ1)-k*IR1
if (idum < 0) idum = idum+IM1
k = idum2/IQ2
idum2 = IA2*(idum2-k*IQ2)-k*IR2
if (idum2 < 0) idum2 = idum2+IM2
j = 1+iy/NDIV
iy = iv(j)-idum2
iv(j) = idum
if (iy < 1) iy = iy+IMM1
ranf = min(AM*iy,RNMX)

return

end function ranf
