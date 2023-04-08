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
        subroutine Exp1 (A,B,dimi,dim2,dime)
!
!       this routine do:
!       expand A(i,pq) -> B(i,p,q)
!       for matrices, where A(i,pq)=A(i,qp)
!
!       parameter description:
!       A       - input matrix (I)
!       B       - outpun matrix (O)
!       dimi    - first dimension in A,B (I)
!       dim2    - # of pq (I)
!       dimf    - # of p (also q) (I)
!
        implicit none
        integer dimi,dim2,dime
        real*8 A(1:dimi,1:dim2)
        real*8 B(1:dimi,1:dime,1:dime)
!
!       help variables
        integer i,p,q,pq
!
        pq=0
        do p=1,dime
        do q=1,p
        pq=pq+1
!
          do i=1,dimi
            B(i,p,q)=A(i,pq)
          end do
!
          do i=1,dimi
            B(i,q,p)=A(i,pq)
          end do
!
        end do
        end do
!
        return
        end
