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
       subroutine defvhlp22 (r1,v,dimr1a,dimr1c,                        &
     & dimvab,dimva,dimvc,adda,addc)
!
!     this routine do
!     V(ab,c)xxx = R1(a,c,b)-R1(b,c,a) x=a,b
!     for syma=symb<symc
!
!     r1        - r1 matrix (I)
!     v        - v matrix (O)
!     dimr1a         - dimension of a (b) in R1 (I)
!     dimr1c         - dimension of c in R1 (I)
!     dimvab        - dimension of ab in V (I)
!     dimva        - dimension of a (b) in V (I)
!     dimvc        - dimension of c in V (I)
!     adda    - additional constat to a (b) (I)
!     addc    - additional constat to c (I)
!
#include "t31.fh"
       integer dimr1a,dimr1c,dimvab,dimva,dimvc,adda,addc
       real*8 r1(1:dimr1a,1:dimr1c,1:dimr1a)
       real*8 v(1:dimvab,1:dimvc)
!
!     help variables
!
       integer a,b,c,ab,ab0,ar1,cr1

!
       do 100 c=1,dimvc
       cr1=c+addc
       do 101 a=2,dimva
       ar1=a+adda
       ab0=nshf(a)
       do 102 b=1,a-1
       v(ab0+b,c)=r1(ar1,cr1,b+adda)
 102    continue
 101    continue
 100    continue
!
       do 200 a=2,dimva
       ar1=a+adda
       ab0=nshf(a)
       do 201 c=1,dimvc
       cr1=c+addc
       do 202 b=1,a-1
       ab=ab0+b
       v(ab,c)=v(ab,c)-r1(b+adda,cr1,ar1)
 202    continue
 201    continue
 200    continue
!
       return
       end
