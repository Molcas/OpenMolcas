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
