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
       subroutine defvhlp82 (r2,v,dimr2b,dimr2a,dimr2c,                 &
     & dimva,dimvb,dimvc,adda,addc)
!
!     this routine do
!     V(a,b,c)aba = - R2(b,a,c)
!     for syma<symc
!
!     r2        - r2 matrix (I)
!     v        - v matrix (O)
!     dimr2b         - dimension of b in R2 (I)
!     dimr2a         - dimension of a in R2 (I)
!     dimr2c         - dimension of c in R2 (I)
!     dimva        - dimension of a in V (I)
!     dimvb        - dimension of b in V (I)
!     dimvc        - dimension of c in V (I)
!     adda    - additional constat to a (I)
!     addc    - additional constat to c (I)
!
       integer dimr2b,dimr2a,dimr2c
       integer dimva,dimvb,dimvc,adda,addc
       real*8 r2(1:dimr2b,1:dimr2c,1:dimr2a)
       real*8 v(1:dimva,1:dimvb,1:dimvc)
!
!     help variables
!
       integer a,b,c,ar2,cr2
!
!
       do 100 a=1,dimva
       ar2=a+adda
       do 101 c=1,dimvc
       cr2=c+addc
       do 102 b=1,dimvb
       v(a,b,c)=-r2(b,cr2,ar2)
 102    continue
 101    continue
 100    continue
!
       return
       end
