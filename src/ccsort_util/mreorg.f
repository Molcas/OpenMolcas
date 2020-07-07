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
c
c     ----------------------------------
c
       subroutine expmpq (wrk,wrksize,
     & syma,typv3,typm,typp,typq,directyes,
     &                    inverseyes)
c
c     this routine realize reorganization to
c     #3 <m,a p q> with given typv3 and typm,p,q  <-  #2 <a,m|pp,qq>
c     #2 is in shape <a,m|pp,qq> for symm,sympp,symqq with types
c     _a  m,pp,qq - 1,5,5
c     #3 <m a p q> may be antisymetrized or not, two parameters (directyes,
c     inverseyes) can be deduced trom typv3 and typm,p,q and syma
c     but for simplicity these are as input parameters
c     this routine allow to use typv3=0 and 2
c
c     syma      - irrep of a
c     typv3     - typ of final #2 (I)
c     typm,p,q  - types of ind. m,p,q (I)
c     directyes - 1 if direct <pqij> integrals are included (I)
c     inverseyes- 1 if inverse <qpij> integrals are included (I)
c
c     foreingh routines used: grc0
c     ccsort_mv0zero
c
c     it also defines new mapd2,mapi2 corresponding to #2
c
#include "wrk.fh"
#include "reorg.fh"
#include "ccsort.fh"
       integer typv3,typp,typq,typm,syma,directyes,inverseyes
c
c     help variables
c
       integer symp,symq,symm,possv3,lenght
       integer ii,iiv2d,iiv2i,possv2d,possv2i
       integer posst
c
c*    get mapd mapi of <m,a|p,q> as _a(m,p q) into mapd3,mapi3
c
       call ccsort_grc0 (3,typv3,typm,typp,typq,0,syma,
     & poss30,posst,mapd3,mapi3)
c
c*    realize reorganization psb
c
       do 100 ii=1,mapd3(0,5)
c
c*    def parameters of #3
       possv3=mapd3(ii,1)
       lenght=mapd3(ii,2)
       symm=mapd3(ii,3)
       symp=mapd3(ii,4)
       symq=mapd3(ii,5)
c
c*    vanish #3
       call ccsort_mv0zero (lenght,lenght,wrk(possv3))
c
       if (directyes.eq.1) then
c
c**   def possition #2 direct (i.e.
       iiv2d=mapi2(symm,symq,1)
       possv2d=mapd2(iiv2d,1)
c
c**   do #3 <m a p q> <- #2 <a m q p> (i.e. direct)
       call  mreorg (wrk,wrksize,
     & symm,symp,symq,typm,typp,typq,
     & 1,3,2,1,5,5,
     & typv3,possv2d,possv3,1.0d0)
c
       end if
c
       if (inverseyes.eq.1) then
c
c**   def possition #2 inverse (i.e. #2 <symq,symp| symi,symj>)
       iiv2i=mapi2(symm,symp,1)
       possv2i=mapd2(iiv2i,1)
c
c**   do #3 <m a q p> <- - #2 <a m p q> (i.e. inverse)
       call  mreorg (wrk,wrksize,
     & symm,symp,symq,typm,typp,typq,
     & 1,2,3,1,5,5,
     & typv3,possv2i,possv3,-1.0d0)
c
       end if
c
 100    continue
c
       return
       end
c
c     -------------------------------------
c
       subroutine mreorg (wrk,wrksize,
     & symp,symq,symr,typp,typq,typr,
     & posspv2,possqv2,possrv2,typpv2,typqv2,typrv2,
     & typv3,possv20,possv30,fact)
c
c     this routine is up level routine for mreorg1 (also more detailed
c     description can be found there).
c     #2 must be of type 0, #3 can be 0, and 2
c     this routine only prepair some constants, required by ireorg1,
c     that can be deduced form input data - dimp,dimqr,dimt-dimv
c
c     symp-r     - symetries of p-r (I)
c     typp-r     - types of indexes p-r in V2 (I)
c     possp-rv2  - possitions of p-r ind. in V2 (I)
c     typp-rv2   - types of indices, corresponding to p-r in V2 (I)
c     typv3      - type of V3 (0,2) (I)
c     possv20,30 - initial possitions of V2 and V3 in wrk (I)
c     fact       - multiplication factors (usually +-1.0d0) (I)
c
#include "wrk.fh"
#include "reorg.fh"
#include "ccsort.fh"
       integer symp,symq,symr,typp,typq,typr
       integer posspv2,possqv2,possrv2,typpv2,typqv2,typrv2
       integer typv3,possv20,possv30
       real*8 fact
c
c     help variables
c
       integer ind(1:4)
       integer nhelp,mhelp,rc,dimp,dimqr
c
c*    define dimensions of V2
c
       call ireorg2 (symp,typpv2,nhelp,rc)
       ind(posspv2)=nhelp
       call ireorg2 (symq,typqv2,nhelp,rc)
       ind(possqv2)=nhelp
       call ireorg2 (symr,typrv2,nhelp,rc)
       ind(possrv2)=nhelp
c
c*    def dimp,dimqr
c
       call ireorg2 (symp,typp,dimp,rc)
c
       call ireorg2 (symq,typq,nhelp,rc)
       call ireorg2 (symr,typr,mhelp,rc)
c
       if ((typv3.eq.2).and.(symq.eq.symr)) then
       dimqr=(nhelp*(nhelp-1))/2
       else
       dimqr=nhelp*mhelp
       end if
c
c*    use mreorg1
c
       call mreorg1 (symp,symq,symr,typp,typq,typr,
     & posspv2,possqv2,possrv2,typpv2,typqv2,typrv2,
     & typv3,wrk(possv20),wrk(possv30),fact,dimp,dimqr,
     & ind(1),ind(2),ind(3))
c
       return
       end
c
c     ---------
c
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
c
