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
       subroutine t3intpck1 (vint,r,dimv1,dimv2,dimv3,dima,dimbc,
     & symq,symr,syms,nob,nvb)
c
c     this routine pack integral block symi,symq,symr,syms
c     R_i(a,bc) = V_i(b,a,c)
c     for symq(b)=syms(c)
c     and write R block onto proper place of opened t3nam file - lunt3
c
c     vint  - integrals for given symetries for given i (I)
c     r     - final R_i matrix (O)
c     dimv1 - 1-st. dimension of V (I)
c     dimv2 - 2-nd. dimension of V (I)
c     dimv3 - 3-rd. dimension of V (I)
c     dima  - dimension of a in R (I)
c     dimbc - dimension of bc in R (I)
c     symq  - symmetry of q (b) (I)
c     symr  - symmetry of r (a) (I)
c     syms  - symmetry of s (c) (I)
c     nob   - number of beta occupied in each irrep (I)
c     nvb   - number of beta virtuals in each irrep (I)
c
       implicit none
#include "reorg.fh"
#include "files_ccsd.fh"
       integer symq,symr,syms
       integer dimv1,dimv2,dimv3,dima,dimbc
       integer nob(1:8)
       integer nvb(1:8)
       real*8 vint(1:dimv1,1:dimv2,1:dimv3)
       real*8 r(1:dima,1:dimbc)
c
c     help variables
c
       integer a,b,c,bc,adda,length,iaddr
c
c*    if there are no beta virtuals - goto write section
       if (nvb(symq)*nvb(symr)*nvb(syms).eq.0) then
       goto 200
       end if
c
c*    calc additional constant for a
       adda=nob(symr)
c
c*    do packing
c
       bc=0
       do 100 b=nob(symq)+1,nob(symq)+nvb(symq)
       do 101 c=nob(syms)+1,b
       bc=bc+1
       do 102 a=1,nvb(symr)
       r(a,bc)=vint(b,a+adda,c)
 102    continue
 101    continue
 100    continue
c
c*    write section
c
 200    length=dima*dimbc
       if (length.gt.0) then
       iaddr=daddr(lunt3)
       call ddafile (lunt3,1,r(1,1),length,iaddr)
       end if

       return
       end
