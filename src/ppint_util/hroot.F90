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

subroutine hroot(x,nn,dpn,pn1,eps)
! Improves the approximate root  x
! In addition we also obtain
! dpn = derivative of h(n) at x
! pn1 = value of h(n-1) at x

use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: x
integer(kind=iwp), intent(in) :: nn
real(kind=wp), intent(out) :: dpn, pn1
real(kind=wp), intent(in) :: eps
integer(kind=iwp) :: iter
real(kind=wp) :: d, dp, p

! iter = 5 sufficient for 8-byte accuracy up to nn = 7
do iter=1,10
  call hrecur(p,dp,pn1,x,nn)
  d = p/dp
  x = x-d
  if (abs(d) <= eps) exit
end do
dpn = dp

return

end subroutine hroot
