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
       subroutine defvhlp62 (r1,v,dimr1a,dimr1b,dimr1c,                 &
     & dimva,dimvb,dimvc,adda)
!
!     this routine do
!     V(a,b,c)abb = R1(a,b,c)
!     for symb<symc
!
!     r1        - r1 matrix (I)
!     v        - v matrix (O)
!     dimr1a         - dimension of a in R1 (I)
!     dimr1b         - dimension of b in R1 (I)
!     dimr1c         - dimension of c in R1 (I)
!     dimva        - dimension of a in V (I)
!     dimvb        - dimension of b in V (I)
!     dimvc        - dimension of c in V (I)
!     adda    - additional constat to a (I)
!
       integer dimr1a,dimr1b,dimr1c
       integer dimva,dimvb,dimvc,adda
       real*8 r1(1:dimr1a,1:dimr1c,1:dimr1b)
       real*8 v(1:dimva,1:dimvb,1:dimvc)
!
!     help variables
!
       integer a,b,c
!
!
       do 100 b=1,dimvb
       do 101 c=1,dimvc
       do 102 a=1,dimva
       v(a,b,c)=r1(a+adda,c,b)
 102    continue
 101    continue
 100    continue
!
       return
       end
