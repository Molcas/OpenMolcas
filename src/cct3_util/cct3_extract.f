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
c     extract files
c     this package contains following files:
c     ext
c     exth1
c     exth2
c     exth3
c     exth4
c     exth5
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       subroutine ext (wrk,wrksize,
     & nind,exttyp,u,v,x,ssu,ssv,ssx,mapda,mapia,ssa,
     & possb0,mapdb,mapib,ssb,rc)
c
c     thsi routine realize extraction
c
c     A(indA) -> B_u(indB) for given u
c
c     nind   - # of indexex in matrix A  (Input)
c     exttyp - type of extraction :  (Input)
c     1 - A(pqrs) -> B_p (qrs); A(pqr) -> B_p(qr) ; A(pq) -> B_p(q)
c     2 - A(pqrs) -> B_q (prs); A(pqr) -> B_q(p,r); A(pq) -> B_q(p)
c     3 - A(pqrs) -> B_r (pqs); A(pqr) -> B_r(pq)
c     4 - A(pqrs) -> B_s (pqr)
c     5 - A(pqrs) -< A_pq(rs) ; A(pqr) -> B_pq(r)
c     6 - A(pqrs) -> B_qr(p,s); A(pqr) -> B_pq(r)
c     7 - A(pqrs) -> B_rs(pq)
c     8 - A(pqrs) -> B_pqr(s)
c     9 - A(pqrs) -> B_qrs(p)
c
c     u      - value of 1.st fix index (I)
c     v      - value of 2.nd fix index (I)
c     x      - value of 3.rd fix index (I)
c     ssu    - symmetry of 1.st fix index (I)
c     ssv    - symmetry of 2.nd fix index (I)
c     ssx    - symmetry of 3.rd fix index (I)
c     mapda  - direct map matrix corresponding to A  (Input)
c     mapia  - inverse map matrix corresponding to A  (Input)
c     ssa    - overall symetry state  of matrix A  (Input)
c     possb0 - initial possition of matrix B in WRK  (Input)
c     mapdb  - direct map matrix corresponding to B  (Output)
c     mapib  - inverse map matrix corresponding to B  (Output)
c     ssb    - overall symetry state  of matrix B  (Output)
c     rc     - return (error) code  (Output)
c
c     Table of extractions
c
c     nind  exttyp          Operation             Implementation
c     4       1     A(p,q,r,s) -> B _p(q,r,s)          Yes
c     A(pq,r,s)  -> B _p(q,r,s)          Yes
c     A(p,qr,s)  -> B _p(qr,s)           NCI
c     A(p,q,rs)  -> B _p(q,rs)           Yes
c     B(pq,rs)   -> B _p(q,rs)           Yes
c
c     4       2     A(p,q,r,s) -> B _q(p,r,s)          Yes
c     A(pq,r,s)  -> B _q(p,r,s)     No (use 4,1)
c     A(p,qr,s)  -> B _q(p,r,s)          NCI
c     A(p,q,rs)  -> B _q(p,rs)           Yes
c     B(pq,rs)   -> B _q(p,rs)      No (use 4,1)
c
c     4       3     A(p,q,r,s) -> B _r(p,q,s)          Yes
c     A(pq,r,s)  -> B _r(pq,s)           Yes
c     A(p,qr,s)  -> B _r(p,q,s)     No (use 4.2 - NCI)
c     A(p,q,rs)  -> B _r(p,q,s)          Yes
c     B(pq,rs)   -> B _r(pq,s)           Yes
c
c     4       4     A(p,q,r,s) -> B _s(p,q,r)          Yes
c     A(pq,r,s)  -> B _s(pq,r)           Yes
c     A(p,qr,s)  -> B _s(p,qr)           NCI
c     A(p,q,rs)  -> B _s(p,q,r)     No (use 4,3)
c     B(pq,rs)   -> B _s(pq,r)      No (use 4,3)
c
c     4       5     A(p,q,r,s) -> B _p_q (r,s)         Yes
c     A(pq,r,s)  -> B _pq (r,s)          NCI
c     A(p,qr,s)  -> B _p_q (r,s)         NCI
c     A(p,q,rs)  -> B _p_q (rs)          NCI
c     A(pq,rs)   -> B _pq (rs)           Yes
c
c     4       6     any case                           NCI
c
c     4       7     A(p,q,r,s) -> B _r_s (p,q)         Yes
c     A(pq,r,s)  -> B _r_s (pq)          NCI
c     A(p,qr,s)  -> B _r_s (p,q)         NCI
c     A(p,q,rs)  -> B _r_s (rs)          Yes
c     A(pq,rs)   -> B _pq (rs)           Yes
c
c     4      8,9    any case                           NCI
c
c     3   any case  any case                           NCI
c
c     2       1     A(p,q)     -> B _p (q)             Yes
c     A(pq)      -> B _p (q)             NCI
c
c     2       2     A(p,q)     -> B _q (p)             Yes
c     A(pq)      -> B _q (p)             NCI
c
c
#include "t31.fh"
#include "wrk.fh"
c
       integer nind,exttyp,ssa,u,v,x,ssu,ssv,ssx,possb0,ssb,rc
       integer mapda(0:512,1:6),mapdb(0:512,1:6)
       integer mapia(1:8,1:8,1:8),mapib(1:8,1:8,1:8)
c
c     help variables
       integer typa,typb,possa,possb,symp,symq,symr,syms
       integer ia,ib,key,posst
       integer nhelp1,nhelp2,jjind,signum
       integer dimp,dimq,dimr,dims
c
c     To fix some warnings
      symr=0
c0.*  some general tests
c
       if (mapda(0,6).eq.2) then
c     RC=2  : nind=4, typA=2 (NCI)
       rc=2
       return
       end if

c
c0.*  def typa and ssB
c
       typa=mapda(0,6)
       if (exttyp.le.4) then
c     one extetnal index
       ssb=mmul(ssa,ssu)
       else if (exttyp.le.7) then
c     two extetnal indexes
       nhelp1=mmul(ssu,ssv)
       ssb=mmul(ssa,nhelp1)
       else
c     3 external indeses
       nhelp1=mmul(ssu,ssv)
       nhelp2=mmul(nhelp1,ssx)
       ssb=mmul(nhelp2,ssa)
       end if
c
c
c
       if (nind.eq.4) then
c
c4    A(pqrs) -> B(stv)
c
       if (exttyp.eq.1) then
c
c4.1  A(pqrs) -> B_p(qrs)
c
c
       if (typa.eq.0) then
c
c4.1.0**** case A(p,q,r,s) -> B_p(q,r,s) ****
c
c4.1.0.*    get mapdb,mapib
       typb=0
       call cct3_grc0 (3,typb,mapda(0,2),mapda(0,3),mapda(0,4),0,ssb,
     & possb0,posst,mapdb,mapib)
c
       do 410 ib=1,mapdb(0,5)
c
c4.1.0.*      def symmetry of all indices
       symq=mapdb(ib,3)
       symr=mapdb(ib,4)
       syms=mapdb(ib,5)
c
c4.1.0.*      find porper A
       ia=mapia(ssu,symq,symr)
c
c4.1.0.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
c
c4.1.0.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
       nhelp1=dimq*dimr*dims
c
c4.1.0.*      realize extraction
       call exth1 (wrk(possa),wrk(possb),dimp,nhelp1,u,1)
c
 410    continue
c
c
       else if (typa.eq.1) then
c
c4.1.1**** case A(pq,r,s) -> B_p(q,r,s) ****
c
       typb=0
       call cct3_grc0 (3,typb,mapda(0,2),mapda(0,3),mapda(0,4),0,ssb,
     & possb0,posst,mapdb,mapib)
c
       do 411 ib=1,mapdb(0,5)
c
c4.1.1.*      def symmetry of all indices
       symq=mapdb(ib,3)
       symr=mapdb(ib,4)
       syms=mapdb(ib,5)
c
c4.1.1.*      def key su>sq - 1 ; su=sq - 2 ; su<sq - 3
       if (ssu.gt.symq) then
       key=1
       else if (ssu.eq.symq) then
       key=2
       else
       key=3
       end if
c
c4.1.1.*      find porper A
       if (key.lt.3) then
       ia=mapia(ssu,symq,symr)
       else
       ia=mapia(symq,ssu,symr)
       end if
c
c4.1.1.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
c
c4.1.1.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
c
c4.1.1.*      realize extraction
       if (key.eq.1) then
       nhelp1=dimq*dimr*dims
       call exth1 (wrk(possa),wrk(possb),dimp,nhelp1,u,1)
       else if (key.eq.2) then
       nhelp1=dimp*(dimp-1)/2
       nhelp2=dimr*dims
       call exth4 (wrk(possa),wrk(possb),dimp,nhelp1,nhelp2,u)
       else
c     key=3
       nhelp1=dimr*dims
       call exth3 (wrk(possa),wrk(possb),dimq,dimp,nhelp1,u,-1)
       end if
c
 411    continue
c
c
       else if (typa.eq.2) then
c
c4.1.2**** case A(p,qr,s) -> B_p(qr,s) ****
c     NCI
c
c
       else if (typa.eq.3) then
c
c4.1.3**** case A(p,q,rs) -> B_p(q,rs) ****
c
       typb=2
       call cct3_grc0 (3,typb,mapda(0,2),mapda(0,3),mapda(0,4),0,ssb,
     & possb0,posst,mapdb,mapib)
c
       do 413 ib=1,mapdb(0,5)
c
c4.1.3.*      def symmetry of all indices
       symq=mapdb(ib,3)
       symr=mapdb(ib,4)
       syms=mapdb(ib,5)
c
c4.1.3.*      find porper A
       ia=mapia(ssu,symq,symr)
c
c4.1.3.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
c
c4.1.3.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
       if (symr.eq.syms) then
       nhelp1=dimq*dimr*(dimr-1)/2
       else
       nhelp1=dimq*dimr*dims
       end if
c
c4.1.3.*      realize extraction
       call exth1 (wrk(possa),wrk(possb),dimp,nhelp1,u,1)
c
 413    continue
c
c
       else if (typa.eq.4) then
c
c4.1.4**** case A(pq,rs) -> B_p(q,rs) ****
c
       typb=2
       call cct3_grc0 (3,typb,mapda(0,2),mapda(0,3),mapda(0,4),0,ssb,
     & possb0,posst,mapdb,mapib)
c
       do 414 ib=1,mapdb(0,5)
c
c4.1.4.*      def symmetry of all indices
       symq=mapdb(ib,3)
       symr=mapdb(ib,4)
       syms=mapdb(ib,5)
c
c4.1.4.*      def key su>sq - 1 ; su=sq - 2 ; su<sq - 3
       if (ssu.gt.symq) then
       key=1
       else if (ssu.eq.symq) then
       key=2
       else
       key=3
       end if
c
c4.1.4.*      find porper A
       if (key.lt.3) then
       ia=mapia(ssu,symq,symr)
       else
       ia=mapia(symq,ssu,symr)
       end if
c
c4.1.4.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
c
c4.1.4.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
c
c4.1.4.*      realize extraction
c
       if (key.eq.1) then
       if (symr.eq.syms) then
       nhelp1=dimq*dimr*(dimr-1)/2
       else
       nhelp1=dimq*dimr*dims
       end if
       call exth1 (wrk(possa),wrk(possb),dimp,nhelp1,u,1)
c
       else if (key.eq.2) then
       nhelp1=dimp*(dimp-1)/2
       if (symr.eq.syms) then
       nhelp2=dimr*(dimr-1)/2
       else
       nhelp2=dimr*dims
       end if
       call exth4 (wrk(possa),wrk(possb),dimp,nhelp1,nhelp2,u)
c
       else
c     key=3
       if (symr.eq.syms) then
       nhelp1=dimr*(dimr-1)/2
       else
       nhelp1=dimr*dims
       end if
       call exth3 (wrk(possa),wrk(possb),dimq,dimp,nhelp1,u,-1)
       end if
c
 414    continue
c
c
       end if
c
c
c
       else if (exttyp.eq.2) then
c
c4.2  A(pqrs) -> B_q(prs)
c
c
       if ((typa.eq.0).or.(typa.eq.3)) then
c
c4.2.0case A(p,q,r,s) -> B_q(p,r,s)
c
       typb=0
       call cct3_grc0 (3,typb,mapda(0,1),mapda(0,3),mapda(0,4),0,ssb,
     & possb0,posst,mapdb,mapib)
c
       do 420 ib=1,mapdb(0,5)
c
c4.2.0.*      def symmetry of all indices
       symp=mapdb(ib,3)
       symr=mapdb(ib,4)
       syms=mapdb(ib,5)
c
c4.2.0.*      find porper A
       ia=mapia(symp,ssu,symr)
c
c4.2.0.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
c
c4.2.0.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
       nhelp1=dimr*dims
c
c4.2.0.*      realize extraction
       call exth3 (wrk(possa),wrk(possb),dimp,dimq,nhelp1,u,1)
c
 420    continue
c
c
       else if (typa.eq.1) then
c
c4.2.1case A(pq,r,s) -> B_q(p,r,s)
c
c     RC=3 : nind=4, typA=1, exttyp=2 (NI - Use exttyp 1)
       rc=3
       return
c
c
       else if (typa.eq.2) then
c
c4.2.2case A(p,qr,s) -> B_q(p,r,s)
c
c     RC=4 : nind=4, typA=2 (NCI)
       rc=4
       return
c
c
       else if (typa.eq.3) then
c
c4.2.3case A(p,q,rs) -> B_q(p,rs)
c
       typb=2
       call cct3_grc0 (3,typb,mapda(0,1),mapda(0,3),mapda(0,4),0,ssb,
     & possb0,posst,mapdb,mapib)
c
       do 423 ib=1,mapdb(0,5)
c
c4.2.3.*      def symmetry of all indices
       symp=mapdb(ib,3)
       symr=mapdb(ib,4)
       syms=mapdb(ib,5)
c
c4.2.3.*      find porper A
       ia=mapia(symp,ssu,symr)
c
c4.2.3.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
c
c4.2.3.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
       if (symr.eq.syms) then
       nhelp1=dimr*(dimr-1)/2
       else
       nhelp1=dimr*dims
       end if
c
c4.2.3.*      realize extraction
       call exth3 (wrk(possa),wrk(possb),dimp,dimq,nhelp1,u,1)
c
 423    continue
c
c
       else if (typa.eq.4) then
c
c4.2.4case A(pq,r,s) -> B_q(p,r,s)
c
c     RC=5 : nind=4, typA=4, exttyp=2 (NI - Use exttyp 1)
       rc=5
       return
c
c
       end if
c
c
c
       else if (exttyp.eq.3) then
c
c4.3  A(pqrs) -> B_r(pqs)
c
c
       if (typa.eq.0) then
c
c4.3.0case A(p,q,r,s) -> B_r(p,q,s)
c
       typb=0
       call cct3_grc0 (3,typb,mapda(0,1),mapda(0,2),mapda(0,4),0,ssb,
     & possb0,posst,mapdb,mapib)
c
       do 430 ib=1,mapdb(0,5)
c
c4.3.0.*      def symmetry of all indices
       symp=mapdb(ib,3)
       symq=mapdb(ib,4)
       syms=mapdb(ib,5)
c
c4.3.0.*      find porper A
       ia=mapia(symp,symq,ssu)
c
c4.3.0.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
c
c4.3.0.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
       nhelp1=dimp*dimq
c
c4.3.0.*      realize extraction
       call exth3 (wrk(possa),wrk(possb),nhelp1,dimr,dims,u,1)
c
 430    continue
c
c
       else if (typa.eq.1) then
c
c4.3.1case A(pq,r,s) -> B_r(pq,s)
c
       typb=1
       call cct3_grc0 (3,typb,mapda(0,1),mapda(0,2),mapda(0,4),0,ssb,
     & possb0,posst,mapdb,mapib)
c
       do 431 ib=1,mapdb(0,5)
c
c4.3.1.*      def symmetry of all indices
       symp=mapdb(ib,3)
       symq=mapdb(ib,4)
       syms=mapdb(ib,5)
c
c4.3.1.*      find porper A
       ia=mapia(symp,symq,ssu)
c
c4.3.1.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
c
c4.3.1.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
       if (symp.eq.symq) then
       nhelp1=dimp*(dimp-1)/2
       else
       nhelp1=dimp*dimq
       end if
c
c4.3.1.*      realize extraction
       call exth3 (wrk(possa),wrk(possb),nhelp1,dimr,dims,u,1)
c
 431    continue
c
c
       else if (typa.eq.2) then
c
c4.3.2case A(p,qr,s) -> B_r(p,q,s)
c
c     RC=6 : nind=4, typA=2 (NCI)
       rc=6
       return
c
c
       else if (typa.eq.3) then
c
c4.3.3case A(p,q,rs) -> B_r(p,q,s)
c
       typb=0
       call cct3_grc0 (3,typb,mapda(0,1),mapda(0,2),mapda(0,4),0,ssb,
     & possb0,posst,mapdb,mapib)
c
       do 433 ib=1,mapdb(0,5)
c
c4.3.3.*      def symmetry of all indices
       symp=mapdb(ib,3)
       symq=mapdb(ib,4)
       syms=mapdb(ib,5)
c
c4.3.3.*      def key su>ss - 1 ; su=ss - 2 ; su<ss - 3
       if (ssu.gt.syms) then
       key=1
       else if (ssu.eq.syms) then
       key=2
       else
       key=3
       end if
c
c4.3.3.*      find porper A
       if (key.lt.3) then
       ia=mapia(symp,symq,ssu)
       else
       ia=mapia(symp,symq,symr)
       end if
c
c4.3.3.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
c
c4.3.3.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
c
c4.3.3.*      realize extraction
       if (key.eq.1) then
       nhelp1=dimp*dimq
       call exth3 (wrk(possa),wrk(possb),nhelp1,dimr,dims,u,1)
       else if (key.eq.2) then
       nhelp1=dimp*dimq
       nhelp2=dimr*(dimr-1)/2
       call exth5 (wrk(possa),wrk(possb),nhelp1,dimr,nhelp2,u)
       else
c     key=3
       nhelp1=dimp*dimq
       call exth3 (wrk(possa),wrk(possb),nhelp1,dims,dimr,u,-1)
       end if
c
 433    continue
c
c
       else if (typa.eq.4) then
c
c4.3.4case A(pq,rs) -> B_r(pq,s)
c
       typb=1
       call cct3_grc0 (3,typb,mapda(0,1),mapda(0,2),mapda(0,4),0,ssb,
     & possb0,posst,mapdb,mapib)
c
       do 434 ib=1,mapdb(0,5)
c
c4.3.4.*      def symmetry of all indices
       symp=mapdb(ib,3)
       symq=mapdb(ib,4)
       syms=mapdb(ib,5)
c
c4.3.4.*      def key su>ss - 1 ; su=ss - 2 ; su<ss - 3
       if (ssu.gt.syms) then
       key=1
       else if (ssu.eq.syms) then
       key=2
       else
       key=3
       end if
c
c4.3.4.*      find porper A
       if (key.lt.3) then
       ia=mapia(symp,symq,ssu)
       else
       ia=mapia(symp,symq,syms)
       end if
c
c4.3.4.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
c
c4.3.4.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
c
c4.3.4.*      realize extraction
       if (key.eq.1) then
       if (symp.eq.symq) then
       nhelp1=dimp*(dimq-1)/2
       else
       nhelp1=dimp*dimq
       end if
       call exth3 (wrk(possa),wrk(possb),nhelp1,dimr,dims,u,1)
       else if (key.eq.2) then
       if (symp.eq.symq) then
       nhelp1=dimp*(dimq-1)/2
       else
       nhelp1=dimp*dimq
       end if
       nhelp2=dimr*(dimr-1)/2
       call exth5 (wrk(possa),wrk(possb),nhelp1,dimr,nhelp2,u)
       else
c     key=3
       if (symp.eq.symq) then
       nhelp1=dimp*(dimq-1)/2
       else
       nhelp1=dimp*dimq
       end if
       nhelp1=nhelp1*dimr
       call exth2 (wrk(possa),wrk(possb),nhelp1,dims,u,-1)
       end if
c
 434    continue
c
c
       end if
c
c
c
       else if (exttyp.eq.4) then
c
c4.4  case A(pqrs) -> B_s(pqr)
c
c
       if (typa.eq.0) then
c
c4.4.0case A(p,q,r,s) -> B_s(p,q,r)
c
       typb=0
       call cct3_grc0 (3,typb,mapda(0,1),mapda(0,2),mapda(0,3),0,ssb,
     & possb0,posst,mapdb,mapib)
c
       do 440 ib=1,mapdb(0,5)
c
c4.4.0.*      def symmetry of all indices
       symp=mapdb(ib,3)
       symq=mapdb(ib,4)
       symr=mapdb(ib,5)
c
c4.4.0.*      find porper A
       ia=mapia(symp,symq,symr)
c
c4.4.0.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
c
c4.4.0.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
       nhelp1=dimp*dimq*dimr
c
c4.4.0.*      realize extraction
       call exth2 (wrk(possa),wrk(possb),nhelp1,dims,u,1)
c
 440    continue
c
c
       else if (typa.eq.1) then
c
c4.4.1case A(pq,r,s) -> B_s(pq,r)
c
       typb=1
       call cct3_grc0 (3,typb,mapda(0,1),mapda(0,2),mapda(0,3),0,ssb,
     & possb0,posst,mapdb,mapib)
c
       do 441 ib=1,mapdb(0,5)
c
c4.4.1.*      def symmetry of all indices
       symp=mapdb(ib,3)
       symq=mapdb(ib,4)
       symr=mapdb(ib,5)
c
c4.4.1.*      find porper A
       ia=mapia(symp,symq,symr)
c
c4.4.1.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
c
c4.4.1.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
c
c4.4.1.*      realize extraction
       if (symp.eq.symq) then
       nhelp1=dimp*(dimq-1)*dimr/2
       else
       nhelp1=dimp*dimq*dimr
       end if
       call exth2 (wrk(possa),wrk(possb),nhelp1,dims,u,1)
c
 441    continue
c
c
       else if (typa.eq.2) then
c
c4.4.2case A(p,qr,s) -> B_s(p,qr)
c
c     RC=7 : nind=4, typA=2  (NCI)
       rc=7
       return
c
c
       else if (typa.eq.3) then
c
c4.4.3case A(p,q,rs) -> B_s(p,q,r)
c
c     RC=8 : nind=4, typA=3, exttyp=4 (NI - Use exttyp 3)
       rc=8
       return
c
c
       else if (typa.eq.4) then
c
c4.4.4case A(pq,rs) -> B_s(pq,r)
c
c     RC=9 : nind=4, typA=4, exttyp=4 (NI - Use exttyp 3)
       rc=9
       return
c
c
       end if
c
       else if (exttyp.eq.5) then
c
c4.5  case A(pqrs) -> B_pq(rs)
c
       if (typa.eq.0) then
c
c4.5.0case A(p,q,r,s) -> B_p_q(r,s)
c
       typb=0
       call cct3_grc0 (2,typb,mapda(0,3),mapda(0,4),0,0,ssb,
     & possb0,posst,mapdb,mapib)
c
       do 450 ib=1,mapdb(0,5)
c
c4.5.0.*      def symmetry of all indices
       symr=mapdb(ib,3)
       syms=mapdb(ib,4)
c
c4.5.0.*      find porper A
       ia=mapia(ssu,ssv,symr)
c
c4.5.0.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
c
c4.5.0.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
       nhelp1=dimp*dimq
       nhelp2=dimr*dims
       jjind=(v-1)*dimp+u
c
c4.5.0.*      realize extraction
       call exth1 (wrk(possa),wrk(possb),nhelp1,nhelp2,jjind,1)
c
 450    continue
c
       else if (typa.eq.4) then
c
c4.5.4case A(pq,rs) -> B_pq(rs)
c
       typb=1
       call cct3_grc0 (2,typb,mapda(0,3),mapda(0,4),0,0,ssb,
     & possb0,posst,mapdb,mapib)
c
       if (ssu.ge.ssv) then
c4.5.4.1      case symp >= symq
c
       do 4541 ib=1,mapdb(0,5)
c
c4.5.4.1.*      def symmetry of all indices
       symr=mapdb(ib,3)
       syms=mapdb(ib,4)
c
c4.5.4.1.*      find porper A
       ia=mapia(ssu,ssv,symr)
c
c4.5.4.1.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
c
c4.5.4.1.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
c
       if (ssu.eq.ssv) then
       nhelp1=dimp*(dimp-1)/2
       else
       nhelp1=dimp*dimq
       end if
c
       if (symr.eq.syms) then
       nhelp2=dimr*(dimr-1)/2
       else
       nhelp2=dimr*dims
       end if
c
c4.5.4.1.*        def joined index and signum
       signum=1
       if (ssu.eq.ssv)then
       if (u.gt.v) then
       jjind=nshf(u)+v
       else if (u.eq.v) then
       jjind=0
       signum=0
       else
       jjind=nshf(v)+u
       signum=-1
       end if
       else
       jjind=(v-1)*dimp+u
       end if
c
c4.5.4.1.*      realize extraction
       call exth1 (wrk(possa),wrk(possb),nhelp1,nhelp2,jjind,
     & signum)
c
 4541   continue
c
       else
c4.5.4.2      case symp < symq
c
       do 4542 ib=1,mapdb(0,5)
c
c4.5.4.2.*      def symmetry of all indices
       symr=mapdb(ib,3)
       syms=mapdb(ib,4)
c
c4.5.4.2.*      find porper A
       ia=mapia(ssv,ssu,symr)
c
c4.5.4.2.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
c
c4.5.4.2.*      def dimensions
       dimq=dimm(mapda(0,1),mapda(ia,3))
       dimp=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
c
       nhelp1=dimp*dimq
c
       if (symr.eq.syms) then
       nhelp2=dimr*(dimr-1)/2
       else
       nhelp2=dimr*dims
       end if
c
c4.5.4.2.*        def joined index and signum
       signum=-1
       jjind=(u-1)*dimp+v
c
c4.5.4.2.*      realize extraction
       call exth1 (wrk(possa),wrk(possb),nhelp1,nhelp2,jjind,
     & signum)
c
 4542   continue
c
       end if
c
       else
c     RC=10 , nindA=4, exttyp=5, typa is nou 0,4 (NCI)
       rc=10
       return
       end if
c
       else if (exttyp.eq.7) then
c
c4.7  case         A(pqrs) -> B_rs(pq)
c
       if (typa.eq.0) then
c
c4.7.0case A(p,q,r,s) -> B_r_s(p,q)
c
       typb=0
       call cct3_grc0 (2,typb,mapda(0,1),mapda(0,2),0,0,ssb,
     & possb0,posst,mapdb,mapib)
c
       do 470 ib=1,mapdb(0,5)
c
c4.7.0.*      def symmetry of all indices
       symp=mapdb(ib,3)
       symq=mapdb(ib,4)
c
c4.7.0.*      find porper A
       ia=mapia(symp,symq,ssu)
c
c4.7.0.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
c
c4.7.0.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
       nhelp1=dimp*dimq
       nhelp2=dimr*dims
       jjind=(v-1)*dimr+u
c
c4.7.0.*      realize extraction
       call exth2 (wrk(possa),wrk(possb),nhelp1,nhelp2,jjind,1)
c
 470    continue
c
       else if (typa.eq.3) then
c
c4.7.3case A(p,q,rs) -> B_rs(p,q)
c
       typb=2
       call cct3_grc0 (2,typb,mapda(0,1),mapda(0,2),0,0,ssb,
     & possb0,posst,mapdb,mapib)
c
       if (ssu.ge.ssv) then
c4.7.3.1      case symr >= syms
c
       do 4731 ib=1,mapdb(0,5)
c
c4.7.3.1.*      def symmetry of all indices
       symp=mapdb(ib,3)
       symq=mapdb(ib,4)
c
c4.7.3.1.*      find porper A
       ia=mapia(symp,symq,ssu)
c
c4.7.3.1.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
c
c4.7.3.1.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
c
       if (ssu.eq.ssv) then
       nhelp1=dimr*(dimr-1)/2
       else
       nhelp1=dimr*dims
       end if
c
       nhelp2=dimp*dimq
c
c4.7.3.1.*        def joined index and signum
       signum=1
       if (ssu.eq.ssv)then
       if (u.gt.v) then
       jjind=nshf(u)+v
       else if (u.eq.v) then
       jjind=0
       signum=0
       else
       jjind=nshf(v)+u
       signum=-1
       end if
       else
       jjind=(v-1)*dimr+u
       end if
c
c4.7.3.1.*      realize extraction
       call exth2 (wrk(possa),wrk(possb),nhelp2,nhelp1,jjind,
     & signum)
c
 4731   continue
c
       else
c4.7.3.2      case symp < symq
c
       do 4732 ib=1,mapdb(0,5)
c
c4.7.3.2.*      def symmetry of all indices
       symp=mapdb(ib,3)
       symq=mapdb(ib,4)
c
c4.7.3.2.*      find porper A
       ia=mapia(symp,symq,ssv)
c
c4.7.3.2.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
c
c4.7.3.2.*      def dimensions
       dimq=dimm(mapda(0,1),mapda(ia,3))
       dimp=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
c
       nhelp1=dimr*dims
c
       nhelp2=dimp*dimq
c
c4.7.3.2.*        def joined index and signum
       signum=-1
       jjind=(u-1)*dimr+v
c
c4.7.3.2.*      realize extraction
       call exth2 (wrk(possa),wrk(possb),nhelp2,nhelp1,jjind,
     & signum)
c
 4732   continue
c
       end if
c
c
       else if (typa.eq.4) then
c
c4.7.4case A(pq,rs) -> B_rs(pq)
c
       typb=1
       call cct3_grc0 (2,typb,mapda(0,1),mapda(0,2),0,0,ssb,
     & possb0,posst,mapdb,mapib)
c
       if (ssu.ge.ssv) then
c4.7.4.1      case symr >= syms
c
       do 4741 ib=1,mapdb(0,5)
c
c4.7.4.1.*      def symmetry of all indices
       symp=mapdb(ib,3)
       symq=mapdb(ib,4)
c
c4.7.4.1.*      find porper A
       ia=mapia(symp,symq,ssu)
c
c4.7.4.1.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
c
c4.7.4.1.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
c
       if (ssu.eq.ssv) then
       nhelp1=dimr*(dimr-1)/2
       else
       nhelp1=dimr*dims
       end if
c
       if (symp.eq.symq) then
       nhelp2=dimp*(dimp-1)/2
       else
       nhelp2=dimp*dimq
       end if
c
c4.7.4.1.*        def joined index and signum
       signum=1
       if (ssu.eq.ssv)then
       if (u.gt.v) then
       jjind=nshf(u)+v
       else if (u.eq.v) then
       jjind=0
       signum=0
       else
       jjind=nshf(v)+u
       signum=-1
       end if
       else
       jjind=(v-1)*dimr+u
       end if
c
c4.7.4.1.*      realize extraction
       call exth2 (wrk(possa),wrk(possb),nhelp2,nhelp1,jjind,
     & signum)
c
 4741   continue
c
       else
c4.7.4.2      case symp < symq
c
       do 4742 ib=1,mapdb(0,5)
c
c4.7.4.2.*      def symmetry of all indices
       symp=mapdb(ib,3)
       symq=mapdb(ib,4)
c
c4.7.4.2.*      find porper A
       ia=mapia(symp,symq,ssv)
c
c4.7.4.2.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
c
c4.7.4.2.*      def dimensions
       dimq=dimm(mapda(0,1),mapda(ia,3))
       dimp=dimm(mapda(0,2),mapda(ia,4))
       dimr=dimm(mapda(0,3),mapda(ia,5))
       dims=dimm(mapda(0,4),mapda(ia,6))
c
       nhelp1=dimr*dims
c
       if (symp.eq.symq) then
       nhelp2=dimp*(dimp-1)/2
       else
       nhelp2=dimp*dimq
       end if
c
c4.7.4.2.*        def joined index and signum
       signum=-1
       jjind=(u-1)*dimr+v
c
c4.7.4.2.*      realize extraction
       call exth2 (wrk(possa),wrk(possb),nhelp2,nhelp1,jjind,
     & signum)
c
 4742   continue
c
       end if
c
       else
c     RC=11 , nindA=4, exttyp=7, typa is nou 0,4 (NCI)
       rc=11
       return
       end if
c
       else
c     RC=12  nindA=4, exttyp 1,2,3,4,5 or 7 (Stup/NCI)
       rc=12
       return
c
       end if
c
       else if (nind.eq.2) then
c
       if (exttyp.eq.1) then
c
c2.1  case A(p q) -> A _p(q)
c
       if (typa.eq.0) then
c
c2.1.0case A(p,q) -> B _p(q)
c
       typb=0
       call cct3_grc0 (1,typb,mapda(0,2),0,0,0,ssb,
     & possb0,posst,mapdb,mapib)
c
       do 210 ib=1,mapdb(0,5)
c
c2.1.0.*      def symmetry of all indices
       symq=mapdb(ib,3)
c
c2.1.0.*      find porper A
       ia=mapia(ssu,1,1)
c
c2.1.0.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
c
c2.1.0.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
c
c2.1.0.*      realize extraction
       call exth1 (wrk(possa),wrk(possb),dimp,dimq,u,1)
c
 210    continue
c
       else
c     RC=13  : nindA=2, exttyp=1, typa is not 0 (NCI)
       rc=13
       return
       end if
c
       else if (exttyp.eq.2) then
c
c2.2  case A(p q) -> A _q(p)
c
       if (typa.eq.0) then
c
c2.2.0case A(p,q) -> B _q(p)
c
       typb=0
       call cct3_grc0 (1,typb,mapda(0,1),0,0,0,ssb,
     & possb0,posst,mapdb,mapib)
c
       do 220 ib=1,mapdb(0,5)
c
c2.2.0.*      def symmetry of all indices
       symp=mapdb(ib,3)
c
c2.2.0.*      find porper A
       ia=mapia(symp,1,1)
c
c2.2.0.*      def possitions of A and B
       possa=mapda(ia,1)
       possb=mapdb(ib,1)
c
c2.2.0.*      def dimensions
       dimp=dimm(mapda(0,1),mapda(ia,3))
       dimq=dimm(mapda(0,2),mapda(ia,4))
c
c2.2.0.*      realize extraction
       call exth2 (wrk(possa),wrk(possb),dimp,dimq,u,1)
c
 220    continue
c
c
       else
c     RC=14  : nindA=2, exttyp=2, typa is not 0 (NCI)
       rc=14
       return
       end if
c
       else
c     RC=15  : nindA=2, exttyp is not 1,2 (Stup)
       rc=15
       return
c
       end if
c
       else
c     RC=16  nindA is not 2 or 4 (NCI)
       rc=16
       return
c
       end if
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(x)
       end
c
c     -------------------------------------------
c
       subroutine exth1 (a,b,dimp,dimq,p,nfact)
c
c     this routine extract A(p,q) -> B_p(q)
c
c     a     - matrxi a (Input)
c     b     - matrix b (Output)
c     dimp  - dimension of p (Input)
c     dimq  - dimension of q (Input)
c     p     - value of index p (Input)
c     nfact - sign (+-1,0) (Input)
c
       integer dimp,dimq,p,nfact
       real*8 a(1:dimp,1:dimq)
       real*8 b(1:dimq)
c
c     help variables
c
       integer q
c
       if (nfact.eq.1) then
       do 10 q=1,dimq
       b(q)=a(p,q)
 10     continue
       else if (nfact.eq.-1) then
       do 20 q=1,dimq
       b(q)=a(p,q)
 20     continue
       else if (nfact.eq.0) then
       do 30 q=1,dimq
       b(q)=0.0d0
 30     continue
       end if
c
       return
       end
c
c     -------------------------------------------
c
       subroutine exth2 (a,b,dimp,dimq,q,nfact)
c
c     this routine extract A(p,q) -> B_q(p)
c
c     a     - matrxi a (Input)
c     b     - matrix b (Output)
c     dimp  - dimension of p (Input)
c     dimq  - dimension of q (Input)
c     q     - value of index q (Input)
c     nfact - sign (+-1,0) (Input)
c
       integer dimp,dimq,q,nfact
       real*8 a(1:dimp,1:dimq)
       real*8 b(1:dimp)
c
c     help variables
c
       integer p
c
       if (nfact.eq.1) then
       do 10 p=1,dimp
       b(p)=a(p,q)
 10     continue
       else if (nfact.eq.-1) then
       do 20 p=1,dimp
       b(p)=-a(p,q)
 20     continue
       else if (nfact.eq.0) then
       do 30 p=1,dimp
       b(p)=0.0d0
 30     continue
       end if
c
       return
       end
c
c     -------------------------------------------
c
       subroutine exth3 (a,b,dimp,dimq,dimr,q,nfact)
c
c     this routine extract A(p,q,r) -> B_q(p,r)
c
c     a     - matrxi a (Input)
c     b     - matrix b (Output)
c     dimp  - dimension of p (Input)
c     dimq  - dimension of q (Input)
c     dimr  - dimension of r (Input)
c     q     - value of index q (Input)
c     nfact - sign (+-1,0) (Input)
c
       integer dimp,dimq,dimr,q,nfact
       real*8 a(1:dimp,1:dimq,1:dimr)
       real*8 b(1:dimp,1:dimr)
c
c     help variables
c
       integer p,r
c
       if (nfact.eq.1) then
       do 10 r=1,dimr
       do 10 p=1,dimp
       b(p,r)=a(p,q,r)
 10     continue
       else if (nfact.eq.-1) then
       do 20 r=1,dimr
       do 20 p=1,dimp
       b(p,r)=-a(p,q,r)
 20     continue
       else if (nfact.eq.0) then
       do 30 r=1,dimr
       do 30 p=1,dimp
       b(p,r)=0.0d0
 30     continue
       end if

c
       return
       end
c
c     -------------------------------------------
c
c
       subroutine exth4 (a,b,dimp,dimpq,dimr,p)
c
c     this routine extract A(pq,r) -> B_p(q,r)
c
c     a     - matrxi a (Input)
c     b     - matrix b (Output)
c     dimp  - dimension of p (and also q) (Input)
c     dimpq - dimension of pq (Input)
c     dimr  - dimension of r (Input)
c     p     - value of index p (Input)
c
#include "t31.fh"
       integer dimp,dimpq,dimr,p
       real*8 a(1:dimpq,1:dimr)
       real*8 b(1:dimp,dimr)
c
c     help variables
c
       integer q,r,qp,pq0
c
       if (p.eq.0) then
       return
       end if
c
c     q>p part
       if (p.gt.1) then
       pq0=nshf(p)
       do 20 r=1,dimr
       do 20 q=1,p-1
       b(q,r)=a(pq0+q,r)
 20     continue
       end if
c
c     q=p part
       do 40 r=1,dimr
       b(p,r)=0.0d0
 40     continue
c
c     r<p part
       if (p.lt.dimp) then
       do 60 r=1,dimr
       do 60 q=p+1,dimp
       qp=nshf(q)+p
       b(q,r)=-a(qp,r)
 60     continue
       end if
c
       return
       end
c
c     -------------------------------------------
c
       subroutine exth5 (a,b,dimp,dimq,dimqr,q)
c
c     this routine extract A(p,qr) -> B_q(p,r)
c
c     a     - matrxi a (Input)
c     b     - matrix b (Output)
c     dimp  - dimension of p (Input)
c     dimq  - dimension of q (and also r) (Input)
c     dimqr - dimension of qr (Input)
c     q     - value of index q (Input)
c
#include "t31.fh"
       integer dimp,dimq,dimqr,q
       real*8 a(1:dimp,1:dimqr)
       real*8 b(1:dimp,dimq)
c
c     help variables
c
       integer r,p,rq,qr0,qr
c
       if (q.eq.0) then
       return
       end if
c
c     r>q part
       if (q.gt.1) then
       qr0=nshf(q)
       do 20 r=1,q-1
       qr=qr0+r
       do 20 p=1,dimp
       b(p,r)=a(p,qr)
 20     continue
       end if
c
c     r=q part
       do 40 p=1,dimp
       b(p,q)=0.0d0
 40     continue
c
c     r<p part
       if (q.lt.dimq) then
       do 60 r=q+1,dimq
       rq=nshf(r)+q
       do 60 p=1,dimp
       b(p,r)=-a(p,rq)
 60     continue
       end if
c
       return
       end
c
c     -------------------------------------------
c
