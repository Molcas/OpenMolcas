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

subroutine Exp2(A,B,dima,dimb,dim_2,dime)
! this routine does:
! expand A(a,b,pq) -> B(a,b,p,q)
! for matrices, where A(a,b,pq)=A(a,b,qp)
!
! parameter description:
! A       - input matrix (I)
! B       - outpun matrix (O)
! dima,b  - first dimension in A,B (I)
! dim_2   - # of pq (I)
! dimf    - # of p (also q) (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: dima, dimb, dim_2, dime
real(kind=wp) :: A(dima,dimb,dim_2), B(dima,dimb,dime,dime)
integer(kind=iwp) :: a1, b1, p, pq, q

pq = 0
do p=1,dime
  do q=1,p
    pq = pq+1

    do b1=1,dimb
      do a1=1,dima
        B(a1,b1,p,q) = A(a1,b1,pq)
      end do
    end do

    do b1=1,dimb
      do a1=1,dima
        B(a1,b1,q,p) = A(a1,b1,pq)
      end do
    end do

  end do
end do

return

end subroutine Exp2
