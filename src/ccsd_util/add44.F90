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

subroutine add44(a,b,r,dimp,dimqr,dimq,fact)
! this routine does:
! B(p,qr) <-- fact * A(p,q) for given r

use ccsd_global, only: nshf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: r, dimp, dimqr, dimq
real(kind=wp), intent(in) :: a(dimp,dimq), fact
real(kind=wp), intent(inout) :: b(dimp,dimqr)
integer(kind=iwp) :: q, qr, rq

if (r /= 1) then

  rq = nshf(r)
  b(:,rq+1:rq+r-1) = b(:,rq+1:rq+r-1)-fact*a(:,1:r-1)

end if

if (r /= dimq) then

  do q=r+1,dimq
    qr = nshf(q)+r
    b(:,qr) = b(:,qr)+fact*a(:,q)
  end do

end if

return

end subroutine add44
