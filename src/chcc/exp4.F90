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

subroutine Exp4(A,B,dima2,dima,dimp2,dimp)
! this routine does:
! expand A(ab,pq) -> B(b,a,p,q)
! for matrices, where A(ab,pq)=A(ab,qp)=A(ba,pq)=A(ba,qp)
!  N.B. kvajto odflaknute
!
! parameter description:
! A       - input matrix (I)
! B       - outpun matrix (O)

implicit none
integer dima2, dima, dimp2, dimp
real*8 A(1:dima2,1:dimp2)
real*8 B(1:dima,1:dima,1:dimp,1:dimp)
! help variables
integer a1, b1, ab, p, q, pq
!
pq = 0
do p=1,dimp
  do q=1,p
    pq = pq+1

    ab = 0
    do a1=1,dima
      do b1=1,a1
        ab = ab+1
        B(a1,b1,p,q) = A(ab,pq)
        B(a1,b1,q,p) = A(ab,pq)
        B(b1,a1,p,q) = A(ab,pq)
        B(b1,a1,q,p) = A(ab,pq)
      end do
    end do

  end do
end do

return

end subroutine Exp4
