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
       subroutine t3intpck2 (vint,r,dimv1,dimv2,dimv3,dima,dimb,dimc,   &
     & symq,symr,syms,nob,nvb)
!
!     this routine pack integral block symi,symq,symr,syms
!     R_i(a,b,c) = V_i(b,a,c)
!     for symq(b)>syms(c)
!     and write R block onto proper place of opened t3nam file - lunt3
!
!     vint  - integrals for given symetries for given i (I)
!     r     - final R_i matrix (O)
!     dimv1 - 1-st. dimension of V (I)
!     dimv2 - 2-nd. dimension of V (I)
!     dimv3 - 3-rd. dimension of V (I)
!     dima  - dimension of a in R (I)
!     dimb  - dimension of b in R (I)
!     dimc  - dimension of c in R (I)
!     symq  - symmetry of q (b) (I)
!     symr  - symmetry of r (a) (I)
!     syms  - symmetry of s (c) (I)
!     nob   - number of beta occupied in each irrep (I)
!     nvb   - number of beta virtuals in each irrep (I)
!
       implicit none
#include "reorg.fh"
#include "files_ccsd.fh"
       integer symq,symr,syms
       integer dimv1,dimv2,dimv3,dima,dimb,dimc
       integer nob(1:8)
       integer nvb(1:8)
       real*8 vint(1:dimv1,1:dimv2,1:dimv3)
       real*8 r(1:dima,1:dimb,1:dimc)
!
!     help variables
!
       integer a,b,c,adda,addb,addc,length,iaddr
!
!*    if there are no beta virtuals - skip
       if (nvb(symq)*nvb(symr)*nvb(syms).eq.0) then
       return
       end if
!
!*    calc additional constants for a,b,c
       adda=nob(symr)
       addb=nob(symq)
       addc=nob(syms)

!
!*    do packing
!
       do 100 c=1,nvb(syms)
       do 101 b=1,nvb(symq)
       do 102 a=1,nvb(symr)
       r(a,b,c)=vint(b+addb,a+adda,c+addc)
 102    continue
 101    continue
 100    continue
!
!*    write section
!
        length=dima*dimb*dimc
       if (length.gt.0) then
       iaddr=daddr(lunt3)
       call ddafile (lunt3,1,r(1,1,1),length,iaddr)
       end if
!
       return
       end
