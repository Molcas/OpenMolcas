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

subroutine add41(a,b,p,dimp,dimpq,dimr,fact)
! this routine does:
! B(pq,r) <-- fact * A(q,r) for given p

use ccsd_global, only: nshf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: p, dimp, dimpq, dimr
real(kind=wp) :: a(dimp,dimr), b(dimpq,dimr), fact
integer(kind=iwp) :: pq, q, qp, r

if (p /= 1) then

  do r=1,dimr
    pq = nshf(p)

    do q=1,p-1
      pq = pq+1
      b(pq,r) = b(pq,r)+fact*a(q,r)
    end do

  end do

end if

if (p /= dimp) then

  do r=1,dimr

    do q=p+1,dimp
      qp = nshf(q)+p
      b(qp,r) = b(qp,r)-fact*a(q,r)
    end do

  end do

end if

return

end subroutine add41
