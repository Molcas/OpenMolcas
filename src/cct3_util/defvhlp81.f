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
       subroutine defvhlp81 (r2,v,dimr2b,dimr2a,dimr2c,
     & dimva,dimvb,dimvc,adda,addc)
c
c     this routine do
c     V(a,b,c)aba = - R2(b,a,c)
c     for syma>symc
c
c     r2        - r2 matrix (I)
c     v        - v matrix (O)
c     dimr2b         - dimension of b in R2 (I)
c     dimr2a         - dimension of a in R2 (I)
c     dimr2c         - dimension of c in R2 (I)
c     dimva        - dimension of a in V (I)
c     dimvb        - dimension of b in V (I)
c     dimvc        - dimension of c in V (I)
c     adda    - additional constat to a (I)
c     addc    - additional constat to c (I)
c
       integer dimr2b,dimr2a,dimr2c
       integer dimva,dimvb,dimvc,adda,addc
       real*8 r2(1:dimr2b,1:dimr2a,1:dimr2c)
       real*8 v(1:dimva,1:dimvb,1:dimvc)
c
c     help variables
c
       integer a,b,c,ar2,cr2
c
c
       do 100 c=1,dimvc
       cr2=c+addc
       do 101 a=1,dimva
       ar2=a+adda
       do 102 b=1,dimvb
       v(a,b,c)=-r2(b,ar2,cr2)
 102    continue
 101    continue
 100    continue
c
       return
       end
