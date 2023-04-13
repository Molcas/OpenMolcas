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

subroutine Exp1(A,B,dimi,dim_2,dime)
! this routine does:
! expand A(i,pq) -> B(i,p,q)
! for matrices, where A(i,pq)=A(i,qp)
!
! parameter description:
! A       - input matrix (I)
! B       - outpun matrix (O)
! dimi    - first dimension in A,B (I)
! dim_2   - # of pq (I)
! dimf    - # of p (also q) (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: dimi, dim_2, dime
real(kind=wp) :: A(dimi,dim_2), B(dimi,dime,dime)
integer(kind=iwp) :: i, p, pq, q

pq = 0
do p=1,dime
  do q=1,p
    pq = pq+1

    do i=1,dimi
      B(i,p,q) = A(i,pq)
    end do

    do i=1,dimi
      B(i,q,p) = A(i,pq)
    end do

  end do
end do

return

end subroutine Exp1
