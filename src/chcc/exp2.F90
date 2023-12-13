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
! B       - output matrix (O)
! dima,b  - first dimension in A,B (I)
! dim_2   - # of pq (I)
! dimf    - # of p (also q) (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimb, dim_2, dime
real(kind=wp), intent(in) :: A(dima,dimb,dim_2)
real(kind=wp), intent(out) :: B(dima,dimb,dime,dime)
integer(kind=iwp) :: p, pq

pq = 0
do p=1,dime

  B(:,:,p,1:p-1) = A(:,:,pq+1:pq+p-1)

  B(:,:,1:p,p) = A(:,:,pq+1:pq+p)

  pq = pq+p
end do

return

end subroutine Exp2
