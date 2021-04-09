************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
        subroutine Exp2 (A,B,dima,dimb,dim2,dime)
c
c       this routine do:
c       expand A(a,b,pq) -> B(a,b,p,q)
c       for matrices, where A(a,b,pq)=A(a,b,qp)
c
c       parameter description:
c       A       - input matrix (I)
c       B       - outpun matrix (O)
c       dima,b  - first dimension in A,B (I)
c       dim2    - # of pq (I)
c       dimf    - # of p (also q) (I)
c
        implicit none
        integer dima,dimb,dim2,dime
        real*8 A(1:dima,1:dimb,1:dim2)
        real*8 B(1:dima,1:dimb,1:dime,1:dime)
c
c       help variables
        integer a1,b1,p,q,pq
c
        pq=0
        do p=1,dime
        do q=1,p
        pq=pq+1
c
          do b1=1,dimb
          do a1=1,dima
            B(a1,b1,p,q)=A(a1,b1,pq)
          end do
          end do
c
          do b1=1,dimb
          do a1=1,dima
            B(a1,b1,q,p)=A(a1,b1,pq)
          end do
          end do
c
        end do
        end do
c
        return
        end
