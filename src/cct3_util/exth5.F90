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

subroutine exth5(a,b,dimp,dimq,dimqr,q)
! this routine extracts A(p,qr) -> B_q(p,r)
!
! a     - matrix a (Input)
! b     - matrix b (Output)
! dimp  - dimension of p (Input)
! dimq  - dimension of q (and also r) (Input)
! dimqr - dimension of qr (Input)
! q     - value of index q (Input)

#include "t31.fh"
integer dimp, dimq, dimqr, q
real*8 a(1:dimp,1:dimqr)
real*8 b(1:dimp,dimq)
! help variables
integer r, p, rq, qr0, qr

if (q == 0) then
  return
end if

! r>q part
if (q > 1) then
  qr0 = nshf(q)
  do r=1,q-1
    qr = qr0+r
    do p=1,dimp
      b(p,r) = a(p,qr)
    end do
  end do
end if

! r=q part
do p=1,dimp
  b(p,q) = 0.0d0
end do

! r<p part
if (q < dimq) then
  do r=q+1,dimq
    rq = nshf(r)+q
    do p=1,dimp
      b(p,r) = -a(p,rq)
    end do
  end do
end if

return

end subroutine exth5
