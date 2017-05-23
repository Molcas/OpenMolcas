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
       subroutine exppqij (wrk,wrksize,
     & typv2,typp,typq,typr,typs,directyes,
     &                     inverseyes)
c
c     this routine realize reorganization to
c     #2 <p q i j> with given typv2 and typp-typs  <-  #1 <pp,qq|i,j>
c     #1 is in shape <pp,qq|i,j> for sympp,symqq,symi>=symj with types
c     pp,qq,i,j -5,5,1,1
c     #2 <p q i j> may be antisymetrized or not, two parameters (directyes,
c     inverseyes) can be deduced trom typv2 and typp-s, but for simplicity
c     thes are as input parameters
c     this routine do not allow to use typv2=2
c
c     typv2     - typ of final #2 (I)
c     typp-s    - types of ind. p-s (I)
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
       integer typv2,typp,typq,typr,typs,directyes,inverseyes
c
c     help variables
c
       integer symp,symq,symi,symj,possv2,lenght
       integer ii,iiv1d,iiv1i,possv1d,possv1i
       integer posst
c
c*    get mapd mapi of <p,q r,s> into mapd2,mapi2
c
       call ccsort_grc0 (4,typv2,typp,typq,typr,typs,1,
     & poss20,posst,mapd2,mapi2)
c
c*    realize reorganization psb
c
       do 100 ii=1,mapd2(0,5)
c
c*    def parameters of #2
       possv2=mapd2(ii,1)
       lenght=mapd2(ii,2)
       symp=mapd2(ii,3)
       symq=mapd2(ii,4)
       symi=mapd2(ii,5)
       symj=mapd2(ii,6)
c
c*    skip this step if lenght=0
       if (lenght.eq.0) then
       goto 100
       end if
c
c*    vanish #2
       call ccsort_mv0zero (lenght,lenght,wrk(possv2))
c
       if (symi.ge.symj) then
c*    case symi>=symj - integrals in #1 are in that shape
c
       if (directyes.eq.1) then
c**   def possition #1 direct (i.e. #1 <symp,symq| symi,symj>)
c     direct integrals are always used
       iiv1d=mapi1(symp,symq,symi)
       possv1d=mapd1(iiv1d,1)
c
c**   do #2 <p q i j> <- #1 <p,q|i,j> (i.e. direct)
c     N.B. Since #1 is  always >= symj
c     so in this case orede of indices in #1 and #2 is the same
       call  ireorg (wrk,wrksize,
     & symp,symq,symi,symj,typp,typq,typr,typs,
     & 1,2,3,4,5,5,1,1,
     & typv2,possv1d,possv2,1.0d0)
c
       end if
c
       if (inverseyes.eq.1) then
c
c**   def possition #1 inverse (i.e. #1 <symq,symp| symi,symj>)
c     inverse integrals are used only if antysymetry is required
       iiv1i=mapi1(symq,symp,symi)
       possv1i=mapd1(iiv1i,1)
c
c**   do #2 <p q i j> <- - #1 <q,p|i,j> (i.e. inverse)
c     N.B. Since #1 is  always >= symj
c     so in this case orede of indices in #1 and #2 are inversed 1<->2
       call  ireorg (wrk,wrksize,
     & symp,symq,symi,symj,typp,typq,typr,typs,
     & 2,1,3,4,5,5,1,1,
     & typv2,possv1i,possv2,-1.0d0)
c
       end if
c
       else
c*    case symi<symj - integrals in #1 are in inverse (symi>=symj) shape
c
       if (directyes.eq.1) then
c
c**   def possition #1 direct (i.e. #1 <symq,symp| symj,symi>)
c     direct integrals are always used
       iiv1d=mapi1(symq,symp,symj)
       possv1d=mapd1(iiv1d,1)
c
c**   do #2 <p q i j> <- #1 <q,p|j,i> (i.e. direct)
c     N.B. Since #1 is  always >= symj
c     so in this case orede of indices in #1 and #2 is inversed 1<->2, 3<->4
       call  ireorg (wrk,wrksize,
     & symp,symq,symi,symj,typp,typq,typr,typs,
     & 2,1,4,3,5,5,1,1,
     & typv2,possv1d,possv2,1.0d0)
c
       end if
c
       if (inverseyes.eq.1) then
c
c**   def possition #1 inverse (i.e. #1 <symp,symq| symj,symi>)
c     inverse integrals are used only if antysymetry is required
       iiv1i=mapi1(symp,symq,symj)
       possv1i=mapd1(iiv1i,1)
c
c**   do #2 <p q i j> <- - #1 <p,q|j,i> (i.e. inverse)
c     N.B. Since #1 is  always >= symj
c     so in this case orede of indices in #1 and #2 are inversed 3<->4
       call  ireorg (wrk,wrksize,
     & symp,symq,symi,symj,typp,typq,typr,typs,
     & 1,2,4,3,5,5,1,1,
     & typv2,possv1i,possv2,-1.0d0)
c
       end if
c
       end if
c
 100    continue
c
       return
       end
c
c     ----------------------------------------
c
       subroutine ireorg (wrk,wrksize,
     & symp,symq,symr,syms,typp,typq,typr,typs,
     & posspv1,possqv1,possrv1,posssv1,typpv1,typqv1,typrv1,typsv1,
     & typv2,possv10,possv20,fact)
c
c     this routine is up level routine for ireorg1 (also more detailed
c     description can be found there).
c     v1 must be of type 0, v2 can be 0,1,3 and 4
c     this routine only prepair some constants, required by ireorg1,
c     that can be deduced form input data - dimpq,dimrs,dimt-dimx
c
c     symp-s     - symetries of p-s (I)
c     typp-s     - types of indexes p-s in V2 (I)
c     possp-sv1  - possitions of p-s ind. in v1 (I)
c     typp-sv1   - types of indices, corresponding to p-s in V1 (I)
c     typv2      - type of V2 (0,1,2,4) (I)
c     possv10,20 - initial possitions of V1 and V2 in wrk (I)
c     fact       - multiplication factors (usually +-1.0d0) (I)
c
#include "wrk.fh"
#include "reorg.fh"
#include "ccsort.fh"
       integer symp,symq,symr,syms,typp,typq,typr,typs
       integer posspv1,possqv1,possrv1,posssv1,typpv1,typqv1,typrv1,
     & typsv1
       integer typv2,possv10,possv20
       real*8 fact
c
c     help variables
c
       integer ind(1:4)
       integer nhelp,mhelp,rc,dimpq,dimrs
c
c*    define dimensions of V1
c
       call ireorg2 (symp,typpv1,nhelp,rc)
       ind(posspv1)=nhelp
       call ireorg2 (symq,typqv1,nhelp,rc)
       ind(possqv1)=nhelp
       call ireorg2 (symr,typrv1,nhelp,rc)
       ind(possrv1)=nhelp
       call ireorg2 (syms,typsv1,nhelp,rc)
       ind(posssv1)=nhelp
c
c*    def dimpq,dimrs
c
       call ireorg2 (symp,typp,nhelp,rc)
       call ireorg2 (symq,typq,mhelp,rc)
c
       if (((typv2.eq.1).or.(typv2.eq.4)).and.(symp.eq.symq)) then
       dimpq=(nhelp*(nhelp-1))/2
       else
       dimpq=nhelp*mhelp
       end if
c
       call ireorg2 (symr,typr,nhelp,rc)
       call ireorg2 (syms,typs,mhelp,rc)
c
       if (((typv2.eq.3).or.(typv2.eq.4)).and.(symr.eq.syms)) then
       dimrs=(nhelp*(nhelp-1))/2
       else
       dimrs=nhelp*mhelp
       end if
c
c*    use ireorg1
c
       call ireorg1 (symp,symq,symr,syms,typp,typq,typr,typs,
     & posspv1,possqv1,possrv1,posssv1,typpv1,typqv1,typrv1,typsv1,
     & typv2,wrk(possv10),wrk(possv20),fact,dimpq,dimrs,
     & ind(1),ind(2),ind(3),ind(4))
c
       return
       end
c
c     ---------
c
       subroutine ireorg1 (symp,symq,symr,syms,typp,typq,typr,typs,
     & posspv1,possqv1,possrv1,posssv1,typpv1,typqv1,typrv1,typsv1,
     & typv2,v1,v2,fact,dimpq,dimrs,dimt,dimu,dimv,dimx)
c
c     this routine map v2(pq,rs) <+ fact . v1 (t,u,v,z)
c     v2 may be of type 0,1,3,4 (typv2) while type v1 is always 0
c     symp-syms and typp-typs are symmetries and types of p-s indexex
c     posspv1-posssv1 are corresponding possitions of p-s indexes in v1
c
c     symp-s     - symetries of p-s (I)
c     typp-s     - types of indexes p-s in V2 (I)
c     possp-sv1  - possitions of p-s ind. in v1 (I)
c     typp-sv1   - types of indices, corresponding to p-s in V1 (I)
c     typv2      - type of V2 (0,1,2,4) (I)
c     v1,v2      - arrays V1 and V2  (I,O)
c     fact       - multiplication factors (usually +-1.0d0) (I)
c     dimpq,rs   - dimensions of V2 (I)
c     dimt-x     - dimensions of V1 (I)
c
c
c     reorg.fh may not be included
#include "ccsort.fh"
       integer symp,symq,symr,syms,typp,typq,typr,typs
       integer posspv1,possqv1,possrv1,posssv1
       integer typpv1,typqv1,typrv1,typsv1,typv2
       integer dimpq,dimrs,dimt,dimu,dimv,dimx
       real*8 v2(1:dimpq,1:dimrs)
       real*8 v1(1:dimt,1:dimu,1:dimv,1:dimx)
       real*8 fact
c
c     help variables
c
       integer p,q,r,s,pq,rs,pup,qup,rup,sup,rc,pqyes,rsyes
       integer paddv1,qaddv1,raddv1,saddv1
       integer ind(1:4)
c
c*    def additive constants
c
       call ireorg3 (symp,typp,typpv1,paddv1,rc)
       call ireorg3 (symq,typq,typqv1,qaddv1,rc)
       call ireorg3 (symr,typr,typrv1,raddv1,rc)
       call ireorg3 (syms,typs,typsv1,saddv1,rc)
c
c*    def sumation limits
c
       call ireorg2 (symp,typp,pup,rc)
       call ireorg2 (symq,typq,qup,rc)
       call ireorg2 (symr,typr,rup,rc)
       call ireorg2 (syms,typs,sup,rc)
c
c*    def pqyes, rsyes (i.e. if there is a reduced sumations)
c
       if ((typv2.eq.1).or.(typv2.eq.4)) then
       if (symp.eq.symq) then
       pqyes=1
       else
       pqyes=0
       end if
       else
       pqyes=0
       end if
c
       if ((typv2.eq.3).or.(typv2.eq.4)) then
       if (symr.eq.syms) then
       rsyes=1
       else
       rsyes=0
       end if
       else
       rsyes=0
       end if
c
c
       if ((pqyes.eq.1).and.(rsyes.eq.1)) then
c
c*    case p>q, r>s
c
       rs=0
       do 100 r=2,rup
       ind(possrv1)=raddv1+r
       do 100 s=1,r-1
       ind(posssv1)=saddv1+s
       rs=rs+1
c
       pq=0
       do 100 p=2,pup
       ind(posspv1)=paddv1+p
       do 100 q=1,p-1
       ind(possqv1)=qaddv1+q
       pq=pq+1
c
       v2(pq,rs)=v2(pq,rs)+fact*v1(ind(1),ind(2),ind(3),ind(4))
c
 100    continue
c
       else if (pqyes.eq.1) then
c
c*    case p>q, r,s
c
       rs=0
       do 200 s=1,sup
       ind(posssv1)=saddv1+s
       do 200 r=1,rup
       ind(possrv1)=raddv1+r
       rs=rs+1
c
       pq=0
       do 200 p=2,pup
       ind(posspv1)=paddv1+p
       do 200 q=1,p-1
       ind(possqv1)=qaddv1+q
       pq=pq+1
c
       v2(pq,rs)=v2(pq,rs)+fact*v1(ind(1),ind(2),ind(3),ind(4))
c
 200    continue
c
       else if (rsyes.eq.1) then
c
c*    case p,q, r>s
c
       rs=0
       do 300 r=2,rup
       ind(possrv1)=raddv1+r
       do 300 s=1,r-1
       ind(posssv1)=saddv1+s
       rs=rs+1
c
       pq=0
       do 300 q=1,qup
       ind(possqv1)=qaddv1+q
       do 300 p=1,pup
       ind(posspv1)=paddv1+p
       pq=pq+1
c
       v2(pq,rs)=v2(pq,rs)+fact*v1(ind(1),ind(2),ind(3),ind(4))
c
 300    continue
c
       else
c
c*    case p,q, r,s
c
       rs=0
       do 400 s=1,sup
       ind(posssv1)=saddv1+s
       do 400 r=1,rup
       ind(possrv1)=raddv1+r
       rs=rs+1
c
       pq=0
       do 400 q=1,qup
       ind(possqv1)=qaddv1+q
       do 400 p=1,pup
       ind(posspv1)=paddv1+p
       pq=pq+1
c
       v2(pq,rs)=v2(pq,rs)+fact*v1(ind(1),ind(2),ind(3),ind(4))
c
 400    continue
c
       end if
c
       return
       end
c
c     ------------
c
       subroutine ireorg3 (symp,typp,typpv1,paddv1,rc)
c
c*    this routine def. constants to be added to index from v2
c     to determine proper index in v1
c     N.B. typp and typpv1 must be compatible (this is testet
c     in this version)
c
c     symp  - irrep of p index (I)
c     typp  - typ of p index in v2 (I)
c     typpv1- typ of corresponding p index in v1 (I)
c     paddv1- constant to be added (O) pv1 = pv2+paddv1
c     rc    - return (error) code (O)
c
#include "ccsort.fh"
       integer symp,typp,typpv1,paddv1,rc
c
       rc=0
c
       if ((typp.eq.1).or.(typp.eq.2)) then
       if ((typpv1.eq.1).or.(typpv1.eq.2).or.(typpv1.eq.5)) then
       paddv1=0
       else
       rc=1
c     RC=1 : typp=1 or 2, incompatible typpv1 (Stup)
       return
       end if
       else if (typp.eq.3) then
       if (typpv1.eq.3) then
       paddv1=0
       else if (typpv1.eq.4) then
       paddv1=nvb(symp)-nva(symp)
       else if (typpv1.eq.5) then
       paddv1=noa(symp)
       else
       rc=2
c     RC=2 : typp=3, incompatible typpv1 (Stup)
       return
       end if
       else if (typp.eq.4) then
       if (typpv1.eq.4) then
       paddv1=0
       else if (typpv1.eq.5) then
       paddv1=nob(symp)
       else
       rc=3
c     RC=3 : typp=4, incompatible typpv1 (Stup)
       return
       end if
       else if (typp.eq.5) then
       if (typpv1.eq.5) then
       paddv1=0
       else
c     RC=4 : typp=5, incompatible typpv1 (Stup)
       return
       end if
       else
       rc=5
c     RC=5 : improper typp (Stup)
       return
       end if
c
       return
       end
c
c     ------------
c
       subroutine ireorg2 (symp,typp,pup,rc)
c
c*    this routine def. sumation limits for given symp and typp
c     i.e. number of indexes for this symmetry and typ
c
c     symp  - irrep of p index (I)
c     typp  - typ of p index in v2 (I)
c     pup   - sumation limit
c     rc    - return (error) code (O)
c
#include "ccsort.fh"
       integer symp,typp,pup,rc
c
       if (typp.eq.1) then
       pup=noa(symp)
       else if (typp.eq.2) then
       pup=nob(symp)
       else if (typp.eq.3) then
       pup=nva(symp)
       else if (typp.eq.4) then
       pup=nvb(symp)
       else if (typp.eq.5) then
       pup=norb(symp)
       else
       rc=1
c     RC=1 : bad typp (Stup)
       return
       end if
c
       return
       end
c
