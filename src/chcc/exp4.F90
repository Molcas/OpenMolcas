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
! B       - output matrix (O)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima2, dima, dimp2, dimp
real(kind=wp), intent(in) :: A(dima2,dimp2)
real(kind=wp), intent(out) :: B(dima,dima,dimp,dimp)
integer(kind=iwp) :: a1, ab, p, pq

pq = 0
do p=1,dimp

  ab = 0
  do a1=1,dima
    B(a1,1:a1-1,p,1:p-1) = A(ab+1:ab+a1-1,pq+1:pq+p-1)
    B(a1,1:a1-1,1:p,p) = A(ab+1:ab+a1-1,pq+1:pq+p)
    B(1:a1,a1,p,1:p-1) = A(ab+1:ab+a1,pq+1:pq+p-1)
    B(1:a1,a1,1:p,p) = A(ab+1:ab+a1,pq+1:pq+p)
    ab = ab+a1
  end do

  pq = pq+p
end do

return

end subroutine Exp4
