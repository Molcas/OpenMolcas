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

subroutine exth4(a,b,dimp,dimpq,dimr,p)
! this routine extracts A(pq,r) -> B_p(q,r)
!
! a     - matrix a (Input)
! b     - matrix b (Output)
! dimp  - dimension of p (and also q) (Input)
! dimpq - dimension of pq (Input)
! dimr  - dimension of r (Input)
! p     - value of index p (Input)

#include "t31.fh"
integer dimp, dimpq, dimr, p
real*8 a(1:dimpq,1:dimr)
real*8 b(1:dimp,dimr)
! help variables
integer q, r, qp, pq0

if (p == 0) then
  return
end if
! q>p part
if (p > 1) then
  pq0 = nshf(p)
  do r=1,dimr
    do q=1,p-1
      b(q,r) = a(pq0+q,r)
    end do
  end do
end if

! q=p part
do r=1,dimr
  b(p,r) = 0.0d0
end do

! r<p part
if (p < dimp) then
  do r=1,dimr
    do q=p+1,dimp
      qp = nshf(q)+p
      b(q,r) = -a(qp,r)
    end do
  end do
end if

return

end subroutine exth4
