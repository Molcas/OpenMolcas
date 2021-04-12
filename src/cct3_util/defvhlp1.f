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
       subroutine defvhlp1 (r1,v,dimr1a,dimr1bc,dimvab,dimvc,add)
c
c     this routine do
c     V(ab,c)xxx = R1(a,bc)-R1(b,ac) x=a,b
c     for syma=symb=symc
c
c     r1        - r1 matrix (I)
c     v        - v matrix (O)
c     dimr1a         - dimension of a (b,c) in R1 (I)
c     dimr1bc        - dimension of bc (ac) in R1 (I)
c     dimvab        - dimension of ab in V (I)
c     dimvc        - dimension of c (a,b) in V (I)
c     add     - additional constat (I)
c     (#of singly occ in syma for alpha, 0 for beta)
c
#include "t31.fh"
       integer dimr1a,dimr1bc,dimvab,dimvc,add
       real*8 r1(1:dimr1a,1:dimr1bc)
       real*8 v(1:dimvab,1:dimvc)
c
c     help variables
c
       integer a,b,c,ab0,ab,ar1,cr1,acr1,bk
c      integer bcr1
c
c      do 100 c=1,dimvc
c      cr1=c+add
c      do 100 a=2,dimvc
c      ar1=a+add
c      ab0=nshf(a)
c      do 100 b=1,a-1
c      bcr1=indab(b+add,cr1)
c      v(ab0+b,c)=r1(ar1,bcr1)
c100    continue
c
       do 101 c=1,dimvc
       cr1=c+add
       do 102 a=2,dimvc
       ar1=a+add
       ab0=nshf(a)
       if (c.le.(a-2)) then
c       b1 <= c1
        bk=cr1*(cr1-1)/2+add
        do b=1,c
c       bcr1=cr1*(cr1-1)/2+add+b
        v(ab0+b,c)=r1(ar1,bk+b)
        end do
c       b1 > c1
        bk=(add+c+1)*(add+c)/2+cr1
        do b=c+1,a-1
c       bcr1=(add+b)*(add+b-1)/2+cr1
        v(ab0+b,c)=r1(ar1,bk)
        bk=bk+add+b
        end do
       else
        bk=cr1*(cr1-1)/2+add
        do b=1,a-1
c       b1 <= c1
c       bcr1=cr1*(cr1-1)/2+add+b
        v(ab0+b,c)=r1(ar1,bk+b)
        end do
       end if
 102   continue
 101   continue
c
       do 200 c=1,dimvc
       cr1=c+add
       do 201 a=2,dimvc
       ab0=nshf(a)
c      acr1=indab(a+add,cr1)
         if ((a+add).gt.cr1) then
           acr1=(a+add)*(a+add-1)/2+cr1
         else
           acr1=cr1*(cr1-1)/2+a+add
         end if
       do 202 b=1,a-1
       ab=ab0+b
       v(ab,c)=v(ab,c)-r1(b+add,acr1)
 202    continue
 201    continue
 200    continue
c
       return
       end
