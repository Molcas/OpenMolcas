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
       subroutine exth1 (a,b,dimp,dimq,p,nfact)
!
!     this routine extract A(p,q) -> B_p(q)
!
!     a     - matrxi a (Input)
!     b     - matrix b (Output)
!     dimp  - dimension of p (Input)
!     dimq  - dimension of q (Input)
!     p     - value of index p (Input)
!     nfact - sign (+-1,0) (Input)
!
       integer dimp,dimq,p,nfact
       real*8 a(1:dimp,1:dimq)
       real*8 b(1:dimq)
!
!     help variables
!
       integer q
!
       if (nfact.eq.1) then
       do 10 q=1,dimq
       b(q)=a(p,q)
 10     continue
       else if (nfact.eq.-1) then
       do 20 q=1,dimq
       b(q)=a(p,q)
 20     continue
       else if (nfact.eq.0) then
       do 30 q=1,dimq
       b(q)=0.0d0
 30     continue
       end if
!
       return
       end
