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

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: r, dimp, dimqr, dimq
real(kind=wp) :: a(dimp,dimq), b(dimp,dimqr), fact
#include "ccsd1.fh"
integer(kind=iwp) :: p, q, qr, rq

if (r /= 1) then

  rq = nshf(r)
  do q=1,r-1
    rq = rq+1

    do p=1,dimp
      b(p,rq) = b(p,rq)-fact*a(p,q)
    end do

  end do

end if

if (r /= dimq) then

  do q=r+1,dimq
    qr = nshf(q)+r
    do p=1,dimp
      b(p,qr) = b(p,qr)+fact*a(p,q)
    end do

  end do

end if

return

end subroutine add44
