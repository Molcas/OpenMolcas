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

subroutine hrecur(pn,dpn,pn1,x,nn)

use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: pn, dpn, pn1
real(kind=wp), intent(in) :: x
integer(kind=iwp), intent(in) :: nn
integer(kind=iwp) :: j
real(kind=wp) :: dp, dp1, dq, fj, fj2, p, p1, q

p1 = One
p = x
dp1 = Zero
dp = One
do j=2,nn
  fj = (j)
  fj2 = Half*(fj-One)
  q = x*p-fj2*p1
  dq = x*dp+p-fj2*dp1
  p1 = p
  p = q
  dp1 = dp
  dp = dq
end do
pn = p
dpn = dp
pn1 = p1

return

end subroutine hrecur
