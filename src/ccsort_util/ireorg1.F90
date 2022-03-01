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
       subroutine ireorg1 (symp,symq,symr,syms,typp,typq,typr,typs,     &
     & posspv1,possqv1,possrv1,posssv1,typpv1,typqv1,typrv1,typsv1,     &
     & typv2,v1,v2,fact,dimpq,dimrs,dimt,dimu,dimv,dimx)
!
!     this routine map v2(pq,rs) <+ fact . v1 (t,u,v,z)
!     v2 may be of type 0,1,3,4 (typv2) while type v1 is always 0
!     symp-syms and typp-typs are symmetries and types of p-s indexex
!     posspv1-posssv1 are corresponding possitions of p-s indexes in v1
!
!     symp-s     - symetries of p-s (I)
!     typp-s     - types of indexes p-s in V2 (I)
!     possp-sv1  - possitions of p-s ind. in v1 (I)
!     typp-sv1   - types of indices, corresponding to p-s in V1 (I)
!     typv2      - type of V2 (0,1,2,4) (I)
!     v1,v2      - arrays V1 and V2  (I,O)
!     fact       - multiplication factors (usually +-1.0d0) (I)
!     dimpq,rs   - dimensions of V2 (I)
!     dimt-x     - dimensions of V1 (I)
!
!
!     reorg.fh may not be included
#include "ccsort.fh"
       integer symp,symq,symr,syms,typp,typq,typr,typs
       integer posspv1,possqv1,possrv1,posssv1
       integer typpv1,typqv1,typrv1,typsv1,typv2
       integer dimpq,dimrs,dimt,dimu,dimv,dimx
       real*8 v2(1:dimpq,1:dimrs)
       real*8 v1(1:dimt,1:dimu,1:dimv,1:dimx)
       real*8 fact
!
!     help variables
!
       integer p,q,r,s,pq,rs,rc,pqyes,rsyes
       integer :: pup=0,qup=0,rup=0,sup=0
       integer :: paddv1=-1,qaddv1=-1,raddv1=-1,saddv1=-1
       integer ind(1:4)
!
!*    def additive constants
!
       call ireorg3 (symp,typp,typpv1,paddv1,rc)
       call ireorg3 (symq,typq,typqv1,qaddv1,rc)
       call ireorg3 (symr,typr,typrv1,raddv1,rc)
       call ireorg3 (syms,typs,typsv1,saddv1,rc)
!
!*    def sumation limits
!
       call ireorg2 (symp,typp,pup,rc)
       call ireorg2 (symq,typq,qup,rc)
       call ireorg2 (symr,typr,rup,rc)
       call ireorg2 (syms,typs,sup,rc)
!
!*    def pqyes, rsyes (i.e. if there is a reduced sumations)
!
       if ((typv2.eq.1).or.(typv2.eq.4)) then
       if (symp.eq.symq) then
       pqyes=1
       else
       pqyes=0
       end if
       else
       pqyes=0
       end if
!
       if ((typv2.eq.3).or.(typv2.eq.4)) then
       if (symr.eq.syms) then
       rsyes=1
       else
       rsyes=0
       end if
       else
       rsyes=0
       end if
!
!
       if ((pqyes.eq.1).and.(rsyes.eq.1)) then
!
!*    case p>q, r>s
!
       rs=0
       do 100 r=2,rup
       ind(possrv1)=raddv1+r
       do 101 s=1,r-1
       ind(posssv1)=saddv1+s
       rs=rs+1
!
       pq=0
       do 102 p=2,pup
       ind(posspv1)=paddv1+p
       do 103 q=1,p-1
       ind(possqv1)=qaddv1+q
       pq=pq+1
!
       v2(pq,rs)=v2(pq,rs)+fact*v1(ind(1),ind(2),ind(3),ind(4))
!
 103    continue
 102    continue
 101    continue
 100    continue
!
       else if (pqyes.eq.1) then
!
!*    case p>q, r,s
!
       rs=0
       do 200 s=1,sup
       ind(posssv1)=saddv1+s
       do 201 r=1,rup
       ind(possrv1)=raddv1+r
       rs=rs+1
!
       pq=0
       do 202 p=2,pup
       ind(posspv1)=paddv1+p
       do 203 q=1,p-1
       ind(possqv1)=qaddv1+q
       pq=pq+1
!
       v2(pq,rs)=v2(pq,rs)+fact*v1(ind(1),ind(2),ind(3),ind(4))
!
 203    continue
 202    continue
 201    continue
 200    continue
!
       else if (rsyes.eq.1) then
!
!*    case p,q, r>s
!
       rs=0
       do 300 r=2,rup
       ind(possrv1)=raddv1+r
       do 301 s=1,r-1
       ind(posssv1)=saddv1+s
       rs=rs+1
!
       pq=0
       do 302 q=1,qup
       ind(possqv1)=qaddv1+q
       do 303 p=1,pup
       ind(posspv1)=paddv1+p
       pq=pq+1
!
       v2(pq,rs)=v2(pq,rs)+fact*v1(ind(1),ind(2),ind(3),ind(4))
!
 303    continue
 302    continue
 301    continue
 300    continue
!
       else
!
!*    case p,q, r,s
!
       rs=0
       do 400 s=1,sup
       ind(posssv1)=saddv1+s
       do 401 r=1,rup
       ind(possrv1)=raddv1+r
       rs=rs+1
!
       pq=0
       do 402 q=1,qup
       ind(possqv1)=qaddv1+q
       do 403 p=1,pup
       ind(posspv1)=paddv1+p
       pq=pq+1
!
       v2(pq,rs)=v2(pq,rs)+fact*v1(ind(1),ind(2),ind(3),ind(4))
!
 403    continue
 402    continue
 401    continue
 400    continue
!
       end if
!
       return
       end
