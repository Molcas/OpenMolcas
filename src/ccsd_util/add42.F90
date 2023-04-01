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

#include "ccsd1.fh"
integer dimq, dimpq, dimr, q
real*8 fact
real*8 b(1:dimpq,1:dimr)
real*8 a(1:dimq,1:dimr)
! help variable
integer pq, qp, r, p

if (q == 1) goto 101

do r=1,dimr
  qp = nshf(q)

  do p=1,q-1
    qp = qp+1
    b(qp,r) = b(qp,r)-fact*a(p,r)
  end do

end do

101 continue
if (q == dimq) then
  return
end if

do r=1,dimr

  do p=q+1,dimq
    pq = nshf(p)+q
    b(pq,r) = b(pq,r)+fact*a(p,r)
  end do

end do

return

end subroutine add42
