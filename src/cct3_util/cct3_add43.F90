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

subroutine cct3_add43(a,b,q,dimp,dimqr,dimr,fact)
! this routine does:
! B(p,qr) <-- fact * A(p,r) for given q

use CCT3_global, only: nshf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: q, dimp, dimqr, dimr
real(kind=wp), intent(in) :: a(dimp,dimr), fact
real(kind=wp), intent(inout) :: b(dimp,dimqr)
integer(kind=iwp) :: qr, r, rq

if (q /= 1) then
  qr = nshf(q)
  b(:,qr+1:qr+q-1) = b(:,qr+1:qr+q-1)+fact*a(:,1:q-1)
end if

do r=q+1,dimr
  rq = nshf(r)+q
  b(:,rq) = b(:,rq)-fact*a(:,r)
end do

return

end subroutine cct3_add43
