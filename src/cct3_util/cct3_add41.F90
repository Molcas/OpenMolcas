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

subroutine cct3_add41(a,b,p,dimp,dimpq,dimr,fact)
! this routine does:
! B(pq,r) <-- fact * A(q,r) for given p

use CCT3_global, only: nshf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: p, dimp, dimpq, dimr
real(kind=wp), intent(in) :: a(dimp,dimr), fact
real(kind=wp), intent(inout) :: b(dimpq,dimr)
integer(kind=iwp) :: pq, q, qp

if (p /= 1) then
  pq = nshf(p)
  b(pq+1:pq+p-1,:) = b(pq+1:pq+p-1,:)+fact*a(1:p-1,:)
end if

do q=p+1,dimp
  qp = nshf(q)+p
  b(qp,:) = b(qp,:)-fact*a(q,:)
end do

return

end subroutine cct3_add41
