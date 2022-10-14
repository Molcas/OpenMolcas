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
       subroutine exth2 (a,b,dimp,dimq,q,nfact)
!
!     this routine extract A(p,q) -> B_q(p)
!
!     a     - matrxi a (Input)
!     b     - matrix b (Output)
!     dimp  - dimension of p (Input)
!     dimq  - dimension of q (Input)
!     q     - value of index q (Input)
!     nfact - sign (+-1,0) (Input)
!
       integer dimp,dimq,q,nfact
       real*8 a(1:dimp,1:dimq)
       real*8 b(1:dimp)
!
!     help variables
!
       integer p
!
       if (nfact.eq.1) then
       do 10 p=1,dimp
       b(p)=a(p,q)
 10     continue
       else if (nfact.eq.-1) then
       do 20 p=1,dimp
       b(p)=-a(p,q)
 20     continue
       else if (nfact.eq.0) then
       do 30 p=1,dimp
       b(p)=0.0d0
 30     continue
       end if
!
       return
       end
