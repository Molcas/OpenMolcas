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
       subroutine defvhlp7 (r1,v,dimr1a,dimr1b,dimr1bc,
     & dimva,dimvb,dimvc,adda)
c
c     this routine do
c     V(a,b,c)abb = R1(a,bc)
c     for symb.eq.symc
c
c     r1        - r1 matrix (I)
c     v        - v matrix (O)
c     dimr1a         - dimension of a in R1 (I)
c     dimr1b         - dimension of b (c) in R1 (I)
c     dimr1bc        - dimension of bc in R1 (I)
c     dimva        - dimension of a in V (I)
c     dimvb        - dimension of b in V (I)
c     dimvc        - dimension of c in V (I)
c     adda    - additional constat to a (I)
c
#include "t31.fh"
       integer dimr1a,dimr1b,dimr1bc
       integer dimva,dimvb,dimvc,adda
       real*8 r1(1:dimr1a,1:dimr1bc)
       real*8 v(1:dimva,1:dimvb,1:dimvc)
c
c     help variables
c
       integer a,b,c,bcr1
c
c
       do 100 c=1,dimvc
       do 101 b=1,dimvb
c      bcr1=indab(b,c)
         if (b.ge.c) then
           bcr1=b*(b-1)/2+c
         else
           bcr1=c*(c-1)/2+b
         end if
       do 102 a=1,dimva
       v(a,b,c)=r1(a+adda,bcr1)
 102    continue
 101    continue
 100    continue
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(dimr1b)
       end
