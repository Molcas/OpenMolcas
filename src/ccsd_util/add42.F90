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

subroutine add42(a,b,q,dimq,dimpq,dimr,fact)
! this routine does:
! B(pq,r) <-- fact * A(p,r) for given q

use ccsd_global, only: nshf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: q, dimq, dimpq, dimr
real(kind=wp), intent(in) :: a(dimq,dimr), fact
real(kind=wp), intent(inout) :: b(dimpq,dimr)
integer(kind=iwp) :: p, pq, qp, r

if (q /= 1) then

  do r=1,dimr
    qp = nshf(q)
    b(qp+1:qp+q-1,r) = b(qp+1:qp+q-1,r)-fact*a(1:q-1,r)
  end do

end if

if (q /= dimq) then

  do r=1,dimr

    do p=q+1,dimq
      pq = nshf(p)+q
      b(pq,r) = b(pq,r)+fact*a(p,r)
    end do

  end do

end if

return

end subroutine add42
