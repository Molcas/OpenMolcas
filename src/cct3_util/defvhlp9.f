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
       subroutine defvhlp9 (r2,v,dimr2b,dimr2a,dimr2ac,                 &
     & dimva,dimvb,dimvc,adda,addc)
!
!     this routine do
!     V(a,b,c)aba = - R2(b,ac)
!     for syma.eq.symc
!
!     r2        - r2 matrix (I)
!     v        - v matrix (O)
!     dimr2b         - dimension of b in R2 (I)
!     dimr2a         - dimension of a (c) in R2 (I)
!     dimr2ac - dimension of ac in R2 (I)
!     dimva        - dimension of a in V (I)
!     dimvb        - dimension of b in V (I)
!     dimvc        - dimension of c in V (I)
!     adda    - additional constat to a (I)
!     addc    - additional constat to c (I)
!
#include "t31.fh"
       integer dimr2b,dimr2a,dimr2ac
       integer dimva,dimvb,dimvc,adda,addc
       real*8 r2(1:dimr2b,1:dimr2ac)
       real*8 v(1:dimva,1:dimvb,1:dimvc)
!
!     help variables
!
       integer a,b,c,cr2,acr2
!
!
       do 100 c=1,dimvc
       cr2=c+addc
       do 101 a=1,dimva
!      acr2=indab(a+adda,cr2)
         if ((a+adda).ge.cr2) then
           acr2=(a+adda)*(a+adda-1)/2+cr2
         else
           acr2=cr2*(cr2-1)/2+a+adda
         end if
       do 102 b=1,dimvb
       v(a,b,c)=-r2(b,acr2)
 102    continue
 101    continue
 100    continue
!
       return
! Avoid unused argument warnings
      if (.false.) call Unused_integer(dimr2a)
       end
