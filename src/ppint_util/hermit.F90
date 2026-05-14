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

subroutine hermit(nn,x,a,eps)
! Calculates the zeros  x(i)  of the nn-th order
! Hermite polynomial. The largest zero will be
! stored in x(1). Also calculates the corresponding
! coefficients  a(i)  of the nn-th order Gauss-Hermite
! quadrature formula of degree 2*nn-1. The factor of
! sqrt(pi) has been removed from the a(i).
! A. H. Stroud & D. Secrest, Gaussian quadrature formulas,
! Prentice-Hall, 1966

use Constants, only: Zero, One, Two, Six, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nn
real(kind=wp), intent(out) :: x(nn), a(nn)
real(kind=wp), intent(in) :: eps
integer(kind=iwp) :: i, n1, n2, ni
real(kind=wp) :: cc, dpn, fn, pn1, s, xt
real(kind=wp), parameter :: sixth = One/Six

fn = real(nn,kind=wp)
n1 = nn-1
n2 = (nn+1)/2
cc = One
s = Zero
do i=1,n1
  s = s+Half
  cc = s*cc
end do
s = (Two*fn+One)**sixth
do i=1,n2
  if (i == 1) then
    ! largest zero
    xt = s**3-1.85575_wp/s
  else if (i == 2) then
    ! second zero
    xt = xt-1.14_wp*fn**0.426_wp/xt
  else if (i == 3) then
    ! third zero
    xt = 1.86_wp*xt-0.86_wp*x(1)
  else if (i == 4) then
    ! fourth zero
    xt = 1.91_wp*xt-0.91_wp*x(2)
  else
    ! all other zeros
    xt = 2.0_wp*xt-x(i-2)
  end if

  call hroot(xt,nn,dpn,pn1,eps)
  x(i) = xt
  a(i) = cc/dpn/pn1
  !write (u6,'(2i4,2es25.17)') nn,i,xt,a(i)
  ni = nn-i+1
  x(ni) = -xt
  a(ni) = a(i)
end do

return

end subroutine hermit
