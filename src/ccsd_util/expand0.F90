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

subroutine expand0(a,b,dimpq,dimp)
! expand a(pq) -> b(p,q)
! assumption: p>q, a(p,q)=-a(q,p)
! RISC version

integer dimpq, dimp
real*8 a(1:dimpq)
real*8 b(1:dimp,1:dimp)
! help variables
integer p, q, pq
real*8 scalar

if (dimp > 1) then

  pq = 0
  do p=2,dimp
    do q=1,p-1
      pq = pq+1
      scalar = a(pq)
      b(p,q) = scalar
      b(q,p) = -scalar
    end do
  end do

end if

do p=1,dimp
  b(p,p) = 0.0d0
end do

return

end subroutine expand0
