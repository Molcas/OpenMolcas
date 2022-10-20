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
integer(kind=iwp) :: q, dimp, dimqr, dimr
real(kind=wp) :: a(dimp,dimr), b(dimp,dimqr), fact
integer(kind=iwp) :: p, qr, r, rq

if (q /= 1) then

  qr = nshf(q)
  do r=1,q-1
    qr = qr+1

    do p=1,dimp
      b(p,qr) = b(p,qr)+fact*a(p,r)
    end do

  end do

end if

if (q /= dimr) then

  do r=q+1,dimr
    rq = nshf(r)+q
    do p=1,dimp
      b(p,rq) = b(p,rq)-fact*a(p,r)
    end do

  end do

end if

return

end subroutine cct3_add43
