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
       subroutine mreorg1 (symp,symq,symr,typp,typq,typr,
     & posspv1,possqv1,possrv1,typpv1,typqv1,typrv1,
     & typv2,v1,v2,fact,dimp,dimqr,dimt,dimu,dimv)
c
c     this routine map v2(p,qr) <+ fact . v1 (t,u,v)
c     v2 may be of type 0,2 (typv2) while type v1 is always 0
c     symp-symr and typp-typr are symmetries and types of p-r indexex
c     posspv1-possrv1 are corresponding possitions of p-r indexes in v1
c     N.B. v1 and v2 have no direct relation to #1 or #2, since here there
c     is no reorg.fh included, v1,v2 corresponds to arbitrary matrices
c
c     symp-r     - symetries of p-r (I)
c     typp-r     - types of indexes p-r in V2 (I)
c     possp-rv1  - possitions of p-r ind. in v1 (I)
c     typp-rv1   - types of indices, corresponding to p-r in V1 (I)
c     typv2      - type of V2 (0,1,2,4) (I)
c     v1,v2      - arrays V1 and V2  (I,O)
c     fact       - multiplication factors (usually +-1.0d0) (I)
c     dimp,qr   - dimensions of V2 (I)
c     dimt-s     - dimensions of V1 (I)
c
c
c     reorg.fh may not be included
#include "ccsort.fh"
       integer symp,symq,symr,typp,typq,typr
       integer posspv1,possqv1,possrv1
       integer typpv1,typqv1,typrv1,typv2
       integer dimp,dimqr,dimt,dimu,dimv
       real*8 v2(1:dimp,1:dimqr)
       real*8 v1(1:dimt,1:dimu,1:dimv)
       real*8 fact
c
c     help variables
c
       integer p,q,r,qr,pup,qup,rup,rc,qryes
       integer paddv1,qaddv1,raddv1
       integer ind(1:4)
c
c*    def additional constants
c
       call ireorg3 (symp,typp,typpv1,paddv1,rc)
       call ireorg3 (symq,typq,typqv1,qaddv1,rc)
       call ireorg3 (symr,typr,typrv1,raddv1,rc)
c
c*    def sumation limits
c
       call ireorg2 (symp,typp,pup,rc)
       call ireorg2 (symq,typq,qup,rc)
       call ireorg2 (symr,typr,rup,rc)
c
c*    def qryes, rsyes (i.e. if there is a reduced sumations)
c
       if (typv2.eq.2) then
       if (symq.eq.symr) then
       qryes=1
       else
       qryes=0
       end if
       else
       qryes=0
       end if
c
c
       if (qryes.eq.1) then
c
c*    case p, q>s
c
       qr=0
       do 100 q=2,qup
       ind(possqv1)=qaddv1+q
       do 101 r=1,q-1
       ind(possrv1)=raddv1+r
       qr=qr+1
c
       do 102 p=1,pup
       ind(posspv1)=paddv1+p
c
       v2(p,qr)=v2(p,qr)+fact*v1(ind(1),ind(2),ind(3))
c
 102    continue
 101    continue
 100    continue
c
       else
c
c*    case p q,r
c
       qr=0
       do 200 r=1,rup
       ind(possrv1)=raddv1+r
       do 201 q=1,qup
       ind(possqv1)=qaddv1+q
       qr=qr+1
c
       do 202 p=1,pup
       ind(posspv1)=paddv1+p
c
       v2(p,qr)=v2(p,qr)+fact*v1(ind(1),ind(2),ind(3))
c
 202    continue
 201    continue
 200    continue
c
       end if
c
       return
       end
