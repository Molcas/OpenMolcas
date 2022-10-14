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
       subroutine defvhlp3 (r1,r2,v,dimr1a,dimr1b,dimr1c,dimr2ac,       &
     & dimva,dimvb,dimvc,adda,addb,addc)
!
!     this routine do
!     V(a,b,c)xxx = R1(a,c,b)-R2(b,ac) x=a,b
!     for syma.ne.symb symc.eq.syma
!
!     r1        - r1 matrix (I)
!     r2        - r2 matrix (I)
!     v        - v matrix (O)
!     dimr1a         - dimension of a in R1 (I)
!     dimr1b         - dimension of b in R1 (I)
!     dimr1c         - dimension of c in R1 (I)
!     dimr2ac - dimension of ac in R2 (I)
!     dimva         - dimension of a in V (I)
!     dimvb        - dimension of b in V (I)
!     dimvc        - dimension of c in V (I)
!     adda    - additional constat to a (I)
!     addb    - additional constat to b (I)
!     addc    - additional constat to c (I)
!
#include "t31.fh"
       integer dimr1a,dimr1b,dimr1c,dimr2ac
       integer dimva,dimvb,dimvc,adda,addb,addc
       real*8 r1(1:dimr1a,1:dimr1c,1:dimr1b)
       real*8 r2(1:dimr1b,1:dimr2ac)
       real*8 v(1:dimva,1:dimvb,1:dimvc)
!
!     help variables
!
       integer a,b,c,br1,cr1,cr2,acr2
!
!
       do 100 b=1,dimvb
       br1=b+addb
       do 101 c=1,dimvc
       cr1=c+addc
       do 102 a=1,dimva
       v(a,b,c)=r1(a+adda,cr1,br1)
 102    continue
 101    continue
 100    continue
!
       do 200 c=1,dimvc
       cr2=c+addc
       do 201 a=1,dimvc
!      acr2=indab(a+adda,cr2)
         if ((a+adda).ge.cr2) then
           acr2=(a+adda)*(a+adda-1)/2+cr2
         else
           acr2=cr2*(cr2-1)/2+a+adda
         end if
       do 202 b=1,dimvb
       v(a,b,c)=v(a,b,c)-r2(b+addb,acr2)
 202    continue
 201    continue
 200    continue
!
       return
       end
