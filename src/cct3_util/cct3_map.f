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
c     this file contains following routines:
c     map
c     cct3_noperm
c     map41
c     map42
c     map31
c     map32
c     map21
c     map22
c     map11
c
c     -----------------------------------------------------------
c
       subroutine cct3_map (wrk,wrksize,
     & nind,p,q,r,s,
     & mapda,mapia,ssa,mapdb,mapib,possb0,posst,rc)
c
c     this routine realize mappings
c
c     B(indb) <-- A(inda)
c     where inda are order of indexes in mtx A and indb = Perm(inda)
c
c     nind   - number of indexes in matrix A (and B)  (Input)
c     p      - possition of 1-st index od mtx A in mtx B  (Input)
c     q      - possition of 2-nd index od mtx A in mtx B  (Input)
c     r      - possition of 3-rd index od mtx A in mtx B  (Input)
c     s      - possition of 4-th index od mtx A in mtx B  (Input)
c     mapda  - direct map matrix corresponding to A  (Input)
c     mapia  - inverse map matrix corresponding to A  (Input)
c     ssa    - overall symetry state  of matrix A  (Input)
c     mapdb  - direct map matrix corresponding to B  (Output)
c     mapib  - inverse map matrix corresponding to B  (Output)
c     possb0 - initial possition of matrix B in WRK  (Input)
c     posst  - final possition of matrix B in WRK (Output, not used yet)
c     rc     - return (error) code  (Output)
c
c
c     The table of implemented permutations
c
c     ninda  typA   p  q  r  s            Operation                 Implemented
c
c     4     0     all 24 comb.    A(1,2,3,4) -> B(p,q,r,s)             Yes
c
c     4     1     1  2  3  4      A(12,3,4)  -> B(12,3,4)              Yes
c     4     1     1  2  4  3      A(12,3,4)  -> B(12,4,3)              Yes
c     4     1     2  3  1  4      A(12,3,4)  -> B(3,12,4)              Yes
c     4     1           4  1      A(12,3,4)  -> B(4,12,3)              Yes
c     4     1     3  4  1  2      A(12,3,4)  -> B(3,4,12)              Yes
c     4     1           2  1      A(12,3,4)  -> B(4,3,12)              Yes
c     4     1     other comb.                                          No
c
c     4     2     1  2  3  4      A(1,23,4)  -> B(1,23,4)              Yes
c     4     2     4  2  3  1      A(1,23,4)  -> B(4,23,1)              Yes
c     4     2     3  1  2  4      A(1,23,4)  -> B(23,1,4)              Yes
c     4     2     4        3      A(1,23,4)  -> B(23,4,1)              Yes
c     4     2     1  3  4  2      A(1,23,4)  -> B(1,4,23)              Yes
c     4     2     2        1      A(1,23,4)  -> B(4,1,23)              Yes
c     4     2     other comb.                                          No
c
c     4     3     1  2  3  4      A(1,2,34)  -> B(1,2,34)              Yes
c     4     3     2  1  3  4      A(1,2,34)  -> B(2,1,34)              Yes
c     4     3     3  4  1  2      A(1,2,34)  -> B(34,1,2)              Yes
c     4     3     4  4            A(1,2,34)  -> B(34,2,1)              Yes
c     4     3     1  4  2  3      A(1,2,34)  -> B(1,34,2)              Yes
c     4     3     4  1            A(1,2,34)  -> B(2,34,1)              Yes
c     4     3     other comb.                                          No
c
c     4     4     1  2  3  4      A(12,34)   -> B(12,34)               Yes
c     4     4     3  4  1  2      A(12,34)   -> B(32,12)               Yes
c     4     4     other comb.                                          No
c
c     3     0     all 6 comb.     A(1,2,3)   -> B(p,q,r)               Yes
c
c     3     1     1  2  3  -      A(12,3)    -> B(12,3)                Yes
c     3     1     2  3  1  -      A(12,3)    -> B(3,12)                Yes
c     3     1     other comb.                                          No
c
c     3     2     1  2  3         A(1,23)    -> B(1,23)                Yes
c     3     2     3  1  2         A(1,23)    -> B(23,1)                Yes
c     3     2     other comb.                                          No
c
c     2     0     all 2 comb.     A(1,2)     -> B(p,q)                 Yes
c
c     2     1     1  2  -  -      A(12)      -> B(12)                  Yes
c     2     1     other comb.                                          No
c
c
#include "t31.fh"
#include "wrk.fh"

c
       integer nind,p,q,r,s,ssa,possb0,posst,rc
       integer mapda(0:512,1:6),mapdb(0:512,1:6)
       integer mapia(1:8,1:8,1:8),mapib(1:8,1:8,1:8)
c
       integer type(1:4),dl(1:4),sa(1:4)
       integer nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,nhelp6,nhelp7
       integer ia,ib,typ,newtyp
c
       rc=0
c
c     some test of correctness
c
       nhelp1=p+q+r+s
       if ((nind.eq.1).and.(nhelp1.eq.1)) then
       else if ((nind.eq.2).and.(nhelp1.eq.3)) then
       else if ((nind.eq.3).and.(nhelp1.eq.6)) then
       else if ((nind.eq.4).and.(nhelp1.eq.10)) then
       else
c     RC=1 - uncorrect indexes
       rc=1
       return
       end if
c
c     ******** No permutation *******
c
       if (nind.eq.1) then
       call cct3_noperm (wrk,wrksize,
     & mapda,mapia,mapdb,mapib,possb0,posst)
       return
       end if
c
       if ((nind.eq.2).and.(p.eq.1).and.(q.eq.2)) then
       call cct3_noperm (wrk,wrksize,
     & mapda,mapia,mapdb,mapib,possb0,posst)
       return
       end if
c
       if ((nind.eq.3).and.(p.eq.1).and.(q.eq.2).and.(r.eq.3)) then
       call cct3_noperm (wrk,wrksize,
     & mapda,mapia,mapdb,mapib,possb0,posst)
       return
       end if
c
       if ((nind.eq.4).and.(p.eq.1).and.(q.eq.2).and.(r.eq.3)
     &     .and.(s.eq.4)) then
       call cct3_noperm (wrk,wrksize,
     & mapda,mapia,mapdb,mapib,possb0,posst)
       return
       end if
c
       if (nind.eq.2) then
c
c     *********** 2 index ***********
c
       typ=mapda(0,6)
c
       if (typ.eq.0) then
c
c2.1  map A(p,q) -> B(q,p)
c
c     get mapdb,mapib
c
       type(p)=mapda(0,1)
       type(q)=mapda(0,2)
       call cct3_grc0 (nind,0,type(1),type(2),0,0,ssa,
     & possb0,posst,mapdb,mapib)
c
       do 210 ia=1,mapda(0,5)
       if (mapda(ia,2).eq.0) goto 210
c
       sa(p)=mapda(ia,3)
       sa(q)=mapda(ia,4)
       ib=mapib(sa(1),1,1)
c
c     possitions of A,B
       nhelp1=mapda(ia,1)
       nhelp2=mapdb(ib,1)
c
c     dimp,dimq
       nhelp3=dimm(mapda(0,1),sa(p))
       nhelp4=dimm(mapda(0,2),sa(q))
       call cct3_map21 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,p,q,1)
c
 210    continue
c
       else
c     RC=2 - bad typ for: nind=2
       rc=2
       end if
c
       else if (nind.eq.3) then
c
c     *********** 3 index ***********
c
       typ=mapda(0,6)
c
       if (typ.eq.0) then
c
c3.1  map A(p,q,r) -> B(p1,q1,r1)
c
c     get mapdb,mapib
c
       type(p)=mapda(0,1)
       type(q)=mapda(0,2)
       type(r)=mapda(0,3)
       call cct3_grc0 (nind,0,type(1),type(2),type(3),0,ssa,
     & possb0,posst,mapdb,mapib)
c
       do 310 ia=1,mapda(0,5)
       if (mapda(ia,2).eq.0) goto 310
c
       sa(p)=mapda(ia,3)
       sa(q)=mapda(ia,4)
       sa(r)=mapda(ia,5)
       ib=mapib(sa(1),sa(2),1)
c
c     possitions of A,B
       nhelp1=mapda(ia,1)
       nhelp2=mapdb(ib,1)
c
c     dimp,dimq,r
       nhelp3=dimm(mapda(0,1),sa(p))
       nhelp4=dimm(mapda(0,2),sa(q))
       nhelp5=dimm(mapda(0,3),sa(r))
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,
     &  nhelp5,p,q,r,1)
c
 310    continue
c
       else if (typ.eq.1) then
c
c3.2  map A(pq,r)
c     => only sophystical order of p,q,r is 2,3,1 i.e. B(r,pq)
c     => newtyp=2
c
       if ((p.eq.2).and.(q.eq.3).and.(r.eq.1)) then
c
c     get mapdb,mapib
c
       type(p)=mapda(0,1)
       type(q)=mapda(0,2)
       type(r)=mapda(0,3)
       call cct3_grc0 (nind,2,type(1),type(2),type(3),0,ssa,
     & possb0,posst,mapdb,mapib)
c
       do 320 ia=1,mapda(0,5)
       if (mapda(ia,2).eq.0) goto 320
c
       sa(p)=mapda(ia,3)
       sa(q)=mapda(ia,4)
       sa(r)=mapda(ia,5)
       ib=mapib(sa(1),sa(2),1)
c
c     possitions of A,B
       nhelp1=mapda(ia,1)
       nhelp2=mapdb(ib,1)
c
c     dimp,dimq,dimr
       nhelp3=dimm(mapda(0,1),sa(p))
       nhelp4=dimm(mapda(0,2),sa(q))
       nhelp5=dimm(mapda(0,3),sa(r))
c
c     dimpq
       if (sa(2).eq.sa(3)) then
       nhelp7=nhelp3*(nhelp3-1)/2
       else
       nhelp7=nhelp3*nhelp4
       end if
       call cct3_map21 (wrk(nhelp1),wrk(nhelp2),nhelp7,nhelp5,2,1,1)
c
 320    continue
c
       else
c     RC=3 - bad order of p,q,r for: nind=3,typ=1
       rc=3
       end if
c
       else if (typ.eq.2) then
c
c3.3  map A(p,qr)
c     => only sophystical order of p,q,r is 3,1,2 i.e. B(qr,p)
c     => typb=1
c
       if ((p.eq.3).and.(q.eq.1).and.(r.eq.2)) then
c
c     get mapdb,mapib
c
       type(p)=mapda(0,1)
       type(q)=mapda(0,2)
       type(r)=mapda(0,3)
       call cct3_grc0 (nind,1,type(1),type(2),type(3),0,ssa,
     & possb0,posst,mapdb,mapib)
c
       do 330 ia=1,mapda(0,5)
       if (mapda(ia,2).eq.0) goto 330
c
       sa(p)=mapda(ia,3)
       sa(q)=mapda(ia,4)
       sa(r)=mapda(ia,5)
       ib=mapib(sa(1),sa(2),1)
c
c     possitions of A,B
       nhelp1=mapda(ia,1)
       nhelp2=mapdb(ib,1)
c
c     dimp,dimq,dimr
       nhelp3=dimm(mapda(0,1),sa(p))
       nhelp4=dimm(mapda(0,2),sa(q))
       nhelp5=dimm(mapda(0,3),sa(r))
c
c     dimpq
       if (sa(1).eq.sa(2)) then
       nhelp7=nhelp4*(nhelp4-1)/2
       else
       nhelp7=nhelp4*nhelp5
       end if
       call cct3_map21 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp7,2,1,1)
c
 330    continue
c
       else
c     RC=4 - bad order of p,q,r for: nind=3,typ=2
       rc=4
       end if
c
       else
c     RC=5 - bad type for: nind=3
       rc=5
       end if
c
       else if (nind.eq.4) then
c
c     *********** 4 index ***********
c
       typ=mapda(0,6)
c
       if (typ.eq.0) then
c
c4.1  map A(p,q,r,s)
c
c     get mapdb,mapib
c
       type(p)=mapda(0,1)
       type(q)=mapda(0,2)
       type(r)=mapda(0,3)
       type(s)=mapda(0,4)
       call cct3_grc0 (nind,0,type(1),type(2),type(3),type(4),ssa,
     & possb0,posst,mapdb,mapib)
c
       do 410 ia=1,mapda(0,5)
       if (mapda(ia,2).eq.0) goto 410
c
       sa(p)=mapda(ia,3)
       sa(q)=mapda(ia,4)
       sa(r)=mapda(ia,5)
       sa(s)=mapda(ia,6)
       ib=mapib(sa(1),sa(2),sa(3))
c
c     possitions of A,B
       nhelp1=mapda(ia,1)
       nhelp2=mapdb(ib,1)
c
c     dimp,dimq,dimr,dims
       nhelp3=dimm(mapda(0,1),sa(p))
       nhelp4=dimm(mapda(0,2),sa(q))
       nhelp5=dimm(mapda(0,3),sa(r))
       nhelp6=dimm(mapda(0,4),sa(s))
       call cct3_map41 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,
     & nhelp5,nhelp6,p,q,r,s,1)
c
 410    continue
c
       else if ((typ.gt.0).and.(typ.lt.4)) then
c
c4.2  map A(p,q,r,s) with one limitation
c
c     def newtype and do corresponding tests
c
       if (typ.eq.1) then
       if ((p.gt.3).or.((q-p).ne.1)) then
c     RC=6 bad order of p,q,r,s, for: nind=4,typ=1
       rc=6
       return
       else
       nhelp1=p
       end if
       else if (typ.eq.2) then
       if ((q.gt.3).or.((r-q).ne.1)) then
c     RC=7 bad order of p,q,r,s, for: nind=4,typ=2
       rc=7
       return
       else
       nhelp1=q
       end if
       else if (typ.eq.3) then
       if ((r.gt.3).or.((s-r).ne.1)) then
c     RC=8 bad order of p,q,r,s, for: nind=4,typ=3
       rc=8
       return
       else
       nhelp1=r
       end if
       end if
       newtyp=nhelp1
c
c     get mapdb,mapib
c
       type(p)=mapda(0,1)
       type(q)=mapda(0,2)
       type(r)=mapda(0,3)
       type(s)=mapda(0,4)
       call cct3_grc0 (nind,newtyp,type(1),type(2),type(3),type(4),ssa,
     & possb0,posst,mapdb,mapib)

       do 420 ia=1,mapda(0,5)
       if (mapda(ia,2).eq.0) goto 420
c
       sa(p)=mapda(ia,3)
       sa(q)=mapda(ia,4)
       sa(r)=mapda(ia,5)
       sa(s)=mapda(ia,6)
       ib=mapib(sa(1),sa(2),sa(3))
c
c     possitions of A,B
       nhelp1=mapda(ia,1)
       nhelp2=mapdb(ib,1)
c
c     dimp,dimq,dimr,dims
       dl(p)=dimm(mapda(0,1),sa(p))
       dl(q)=dimm(mapda(0,2),sa(q))
       dl(r)=dimm(mapda(0,3),sa(r))
       dl(s)=dimm(mapda(0,4),sa(s))
c
       if ((newtyp.eq.1).and.(sa(1).eq.sa(2))) then
c     B(p1q1,r1,s1) case => oldtyp can be 2,3
       nhelp7=dl(1)*(dl(1)-1)/2
c
       if (typ.eq.1) then
c     A(pq,r,s) case
c     A(pq,r,s) -> A(pq,s,r)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),nhelp7,dl(r),
     &  dl(s),1,3,2,1)
       else if (typ.eq.2) then
c     A(p,qr,s) case
       if (p.eq.3) then
c     A (p,qr,s) -> B(qr,p,s)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),dl(p),nhelp7,dl(s),
     &  2,1,3,1)
       else
c     A (p,qr,s) -> B(qr,s,p)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),dl(p),nhelp7,dl(s),
     &  3,1,2,1)
       end if
       else
c     A(p,q,rs) case
       if (p.eq.3) then
c     A (p,q,rs) -> B(rs,p,q)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),dl(p),dl(q),
     &  nhelp7,2,3,1,1)
       else
c     A (p,q,rs) -> B(rs,q,p)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),dl(p),dl(q),
     &  nhelp7,3,2,1,1)
       end if
       end if
c
       else if ((newtyp.eq.2).and.(sa(2).eq.sa(3))) then
c     B(p1,q1r1,s1) case => oldtyp can be 1,3
       nhelp7=dl(2)*(dl(2)-1)/2
       if (typ.eq.1) then
c     A(pq,r,s) case
       if (r.eq.1) then
c     A (pq,r,s) -> B(r,pq,s)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),nhelp7,dl(r),dl(s),
     &  2,1,3,1)
       else
c     A (pq,r,s) -> B(s,pq,p)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),nhelp7,dl(r),dl(s),
     &  2,3,1,1)
       end if
       else if (typ.eq.2) then
c     A(p,qr,s) case
c     A(p,qr,s) -> B(s,qr,p)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),dl(p),nhelp7,dl(s),
     & 3,2,1,1)
       else
c     A(p,q,rs) case
       if (p.eq.1) then
c     A (p,q,rs) -> B(p,rs,q)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),dl(p),dl(q),nhelp7,
     & 1,3,2,1)
       else
c     A (p,q,rs) -> B(q,rs,p)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),dl(p),dl(q),nhelp7,
     & 3,1,2,1)
       end if
       end if
c
       else if ((newtyp.eq.3).and.(sa(3).eq.sa(4))) then
c     B(p1,q1,r1s1) case => oldtyp can be 1,2
       nhelp7=dl(3)*(dl(3)-1)/2
       if (typ.eq.1) then
c     A(pq,r,s) case
       if (r.eq.2) then
c     A (pq,r,s) -> B(r,s,pq)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),nhelp7,dl(r),dl(s),
     & 3,1,2,1)
       else
c     A (pq,r,s) -> B(s,r,pq)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),nhelp7,dl(r),dl(s),
     & 3,2,1,1)
       end if
       else if (typ.eq.2) then
c     A(p,qr,s) case
       if (p.eq.1) then
c     A (p,qr,s) -> B(p,s,qr)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),dl(p),nhelp7,dl(s),
     & 1,3,2,1)
       else
c     A (p,qr,s) -> B(s,p,qr)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),dl(p),nhelp7,dl(s),
     & 2,3,1,1)
       end if
       else
c     A(p,q,rs) case
c     A(p,q,rs) -> A(q,p,rs)
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),dl(p),dl(q),nhelp7,
     &  2,1,3,1)
       end if
c
       else
c     B(p1,q1,r1,s1) case
       call cct3_map41 (wrk(nhelp1),wrk(nhelp2),dl(p),dl(q),dl(r),
     & dl(s),p,q,r,s,1)
       end if
c
 420    continue
c
c
       else if (typ.eq.4) then
c
c4.2  map A(pq,rs) -> B(rs,pq)
c
c     tests
c
       if ((p.ne.3).and.(q.ne.4).and.(r.ne.1).and.(s.ne.2)) then
c     RC=9 : nind=4, typA=4 (bad order of p,q,r,s, Stup)
       rc=9
       return
       end if
c
c     get mapdb,mapib
c
       type(p)=mapda(0,1)
       type(q)=mapda(0,2)
       type(r)=mapda(0,3)
       type(s)=mapda(0,4)
       call cct3_grc0 (nind,4,type(1),type(2),type(3),type(4),ssa,
     & possb0,posst,mapdb,mapib)
c
       do 430 ia=1,mapda(0,5)
       if (mapda(ia,2).eq.0) goto 430
c
       sa(p)=mapda(ia,3)
       sa(q)=mapda(ia,4)
       sa(r)=mapda(ia,5)
       sa(s)=mapda(ia,6)
       ib=mapib(sa(1),sa(2),sa(3))
c
c     possitions of A,B
       nhelp1=mapda(ia,1)
       nhelp2=mapdb(ib,1)
c
c     dimp,dimq,dimr,dims
       nhelp3=dimm(mapda(0,1),sa(p))
       nhelp4=dimm(mapda(0,2),sa(q))
       nhelp5=dimm(mapda(0,3),sa(r))
       nhelp6=dimm(mapda(0,4),sa(s))
c
       if ((sa(p).eq.sa(q)).and.(sa(r).eq.sa(s))) then
c     A(pq,rs) -> B(rs,pq)
       nhelp3=nhelp3*(nhelp3-1)/2
       nhelp5=nhelp5*(nhelp5-1)/2
       call cct3_map21 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp5,2,1,1)
c
       else if (sa(p).eq.sa(q)) then
c     A(pq,r,s) -> B(r,s,pq)
       nhelp3=nhelp3*(nhelp3-1)/2
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp5,nhelp6,
     & 3,1,2,1)
c
       else if (sa(r).eq.sa(s)) then
c     A(p,q,rs) -> B(rs,p,q)
       nhelp5=nhelp5*(nhelp5-1)/2
       call cct3_map31 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,
     & 2,3,1,1)
c
       else
c     A(p,q,r,s) -> B(r,s,p,q)
       call cct3_map41 (wrk(nhelp1),wrk(nhelp2),nhelp3,nhelp4,nhelp5,
     & nhelp6,3,4,1,2,1)
       end if
c
 430    continue
c
       else
c     RC=10 bad type for: nind=4
       rc=10
       end if
c
       else
c
c     *********  more than 4 indexes
c
c     RC=11 bad number of indexes
       rc=11
       end if
c
c
       return
       end
c
c     -----------
c
       subroutine cct3_noperm (wrk,wrksize,
     & mapda,mapia,mapdb,mapib,poss0,posst)
c
c     realize mapping without permutation
c     define mapd,mapi
c
#include "t31.fh"
#include "wrk.fh"
c
       integer poss0,posst
       integer mapda(0:512,1:6),mapdb(0:512,1:6)
       integer mapia(1:8,1:8,1:8),mapib(1:8,1:8,1:8)
c
c     help variables
c
       integer ib,nhelp,i,j,k
c
c     def mapib
c
       do 10 k=1,nsym
       do 10 j=1,nsym
       do 10 i=1,nsym
       mapib(i,j,k)=mapia(i,j,k)
 10     continue
c
c     def initial values
c
       do nhelp=1,6
       mapdb(0,nhelp)=mapda(0,nhelp)
       end do
c
       posst=poss0
       do 100 ib=1,mapda(0,5)
       do nhelp=2,6
       mapdb(ib,nhelp)=mapda(ib,nhelp)
       end do
       mapdb(ib,1)=posst
       posst=posst+mapdb(ib,2)
c
       call cct3_map11 (wrk(mapda(ib,1)),wrk(mapdb(ib,1)),mapda(ib,2),1)
c
 100    continue
c
       return
       end
c
c     --------------------------
c
       subroutine cct3_map41 (a,b,dimp,dimq,dimr,dims,p,q,r,s,nfact)
c
c     maping A(p1,q1,r1,s1) -> nfact* B(p2,q2,r2,s2)
c
       real*8 a(*)
       real*8 b(*)
       integer dim(4)
       integer dimp,dimq,dimr,dims,p,q,r,s,nfact
       dim(p)=dimp
       dim(q)=dimq
       dim(r)=dimr
       dim(s)=dims
       call cct3_map42 (a,b,dimp,dimq,dimr,dims,dim(1),dim(2),dim(3),
     & dim(4),p,q,r,s,nfact)
       return
       end
c
       subroutine cct3_map42 (a,b,dimp,dimq,dimr,dims,dim1,dim2,
     & dim3,dim4,p,q,r,s,nfact)
c
       integer dimp,dimq,dimr,dims,dim1,dim2,dim3,dim4,p,q,r,s,nfact
       real*8 a(1:dimp,1:dimq,1:dimr,1:dims)
       real*8 b(1:dim1,1:dim2,1:dim3,1:dim4)
c     integer index(1:4)
c
       integer pp,qq,rr,ss
c
       if (nfact.eq.1) then
c
c     factor + 1
c
c     do 100 ss=1,dims
c     index(s)=ss
c     do 100 rr=1,dimr
c     index(r)=rr
c     do 100 qq=1,dimq
c     index(q)=qq
c     do 100 pp=1,dimp
c     index(p)=pp
c     b(index(1),index(2),index(3),index(4))=a(pp,qq,rr,ss)
c100  continue
c
       if (p.eq.1) then
c     1***
       if (q.eq.2) then
c     12**
       if (r.eq.3) then
c     123* (4)
       do 111 ss=1,dims
       do 111 rr=1,dimr
       do 111 qq=1,dimq
       do 111 pp=1,dimp
       b(pp,qq,rr,ss)=a(pp,qq,rr,ss)
 111    continue
       else
c     124* (3)
       do 112 ss=1,dims
       do 112 rr=1,dimr
       do 112 qq=1,dimq
       do 112 pp=1,dimp
       b(pp,qq,ss,rr)=a(pp,qq,rr,ss)
 112    continue
       end if
       else if (q.eq.3) then
c     13**
       if (r.eq.2) then
c     132* (4)
       do 113 ss=1,dims
       do 113 rr=1,dimr
       do 113 qq=1,dimq
       do 113 pp=1,dimp
       b(pp,rr,qq,ss)=a(pp,qq,rr,ss)
 113    continue
       else
c     134* (2)
       do 114 ss=1,dims
       do 114 rr=1,dimr
       do 114 qq=1,dimq
       do 114 pp=1,dimp
       b(pp,ss,qq,rr)=a(pp,qq,rr,ss)
 114    continue
       end if
       else if (q.eq.4) then
c     14**
       if (r.eq.2) then
c     142* (3)
       do 115 ss=1,dims
       do 115 rr=1,dimr
       do 115 qq=1,dimq
       do 115 pp=1,dimp
       b(pp,rr,ss,qq)=a(pp,qq,rr,ss)
 115    continue
       else
c     143* (2)
       do 116 ss=1,dims
       do 116 rr=1,dimr
       do 116 qq=1,dimq
       do 116 pp=1,dimp
       b(pp,ss,rr,qq)=a(pp,qq,rr,ss)
 116    continue
       end if
       end if
c
       else if (p.eq.2) then
c     2***
       if (q.eq.1) then
c     21**
       if (r.eq.3) then
c     213* (4)
       do 121 ss=1,dims
       do 121 rr=1,dimr
       do 121 qq=1,dimq
       do 121 pp=1,dimp
       b(qq,pp,rr,ss)=a(pp,qq,rr,ss)
 121    continue
       else
c     214* (3)
       do 122 ss=1,dims
       do 122 rr=1,dimr
       do 122 qq=1,dimq
       do 122 pp=1,dimp
       b(qq,pp,ss,rr)=a(pp,qq,rr,ss)
 122    continue
       end if
       else if (q.eq.3) then
c     23**
       if (r.eq.1) then
c     231* (4)
       do 123 ss=1,dims
       do 123 rr=1,dimr
       do 123 qq=1,dimq
       do 123 pp=1,dimp
       b(rr,pp,qq,ss)=a(pp,qq,rr,ss)
 123    continue
       else
c     234* (1)
       do 124 ss=1,dims
       do 124 rr=1,dimr
       do 124 qq=1,dimq
       do 124 pp=1,dimp
       b(ss,pp,qq,rr)=a(pp,qq,rr,ss)
 124    continue
       end if
       else if (q.eq.4) then
c     24**
       if (r.eq.3) then
c     243* (1)
       do 125 ss=1,dims
       do 125 rr=1,dimr
       do 125 qq=1,dimq
       do 125 pp=1,dimp
       b(ss,pp,rr,qq)=a(pp,qq,rr,ss)
 125    continue
       else
c     241* (3)
       do 126 ss=1,dims
       do 126 rr=1,dimr
       do 126 qq=1,dimq
       do 126 pp=1,dimp
       b(rr,pp,ss,qq)=a(pp,qq,rr,ss)
 126    continue
       end if
       end if
c
       else if (p.eq.3) then
c     3***
       if (q.eq.1) then
c     31**
       if (r.eq.2) then
c     312* (4)
       do 131 ss=1,dims
       do 131 rr=1,dimr
       do 131 qq=1,dimq
       do 131 pp=1,dimp
       b(qq,rr,pp,ss)=a(pp,qq,rr,ss)
 131    continue
       else
c     314* (2)
       do 132 ss=1,dims
       do 132 rr=1,dimr
       do 132 qq=1,dimq
       do 132 pp=1,dimp
       b(qq,ss,pp,rr)=a(pp,qq,rr,ss)
 132    continue
       end if
       else if (q.eq.2) then
c     32**
       if (r.eq.1) then
c     321* (4)
       do 133 ss=1,dims
       do 133 rr=1,dimr
       do 133 qq=1,dimq
       do 133 pp=1,dimp
       b(rr,qq,pp,ss)=a(pp,qq,rr,ss)
 133    continue
       else
c     324* (1)
       do 134 ss=1,dims
       do 134 rr=1,dimr
       do 134 qq=1,dimq
       do 134 pp=1,dimp
       b(ss,qq,pp,rr)=a(pp,qq,rr,ss)
 134    continue
       end if
       else if (q.eq.4) then
c     34**
       if (r.eq.1) then
c     341* (2)
       do 135 ss=1,dims
       do 135 rr=1,dimr
       do 135 qq=1,dimq
       do 135 pp=1,dimp
       b(rr,ss,pp,qq)=a(pp,qq,rr,ss)
 135    continue
       else
c     342* (1)
       do 136 ss=1,dims
       do 136 rr=1,dimr
       do 136 qq=1,dimq
       do 136 pp=1,dimp
       b(ss,rr,pp,qq)=a(pp,qq,rr,ss)
 136    continue
       end if
       end if
c
       else if (p.eq.4) then
c     4***
       if (q.eq.1) then
c     41**
       if (r.eq.3) then
c     413* (2)
       do 141 ss=1,dims
       do 141 rr=1,dimr
       do 141 qq=1,dimq
       do 141 pp=1,dimp
       b(qq,ss,rr,pp)=a(pp,qq,rr,ss)
 141    continue
       else
c     412* (3)
       do 142 ss=1,dims
       do 142 rr=1,dimr
       do 142 qq=1,dimq
       do 142 pp=1,dimp
       b(qq,rr,ss,pp)=a(pp,qq,rr,ss)
 142    continue
       end if
       else if (q.eq.2) then
c     42**
       if (r.eq.1) then
c     421* (3)
       do 143 ss=1,dims
       do 143 rr=1,dimr
       do 143 qq=1,dimq
       do 143 pp=1,dimp
       b(rr,qq,ss,pp)=a(pp,qq,rr,ss)
 143    continue
       else
c     423* (1)
       do 144 ss=1,dims
       do 144 rr=1,dimr
       do 144 qq=1,dimq
       do 144 pp=1,dimp
       b(ss,qq,rr,pp)=a(pp,qq,rr,ss)
 144    continue
       end if
       else if (q.eq.3) then
c     43**
       if (r.eq.1) then
c     431* (2)
       do 145 ss=1,dims
       do 145 rr=1,dimr
       do 145 qq=1,dimq
       do 145 pp=1,dimp
       b(rr,ss,qq,pp)=a(pp,qq,rr,ss)
 145    continue
       else
c     432* (1)
       do 146 ss=1,dims
       do 146 rr=1,dimr
       do 146 qq=1,dimq
       do 146 pp=1,dimp
       b(ss,rr,qq,pp)=a(pp,qq,rr,ss)
 146    continue
       end if
       end if
c
       end if
c
       else
c
c     factor = -1
c
c     do 200 ss=1,dims
c     index(s)=ss
c     do 200 rr=1,dimr
c     index(r)=rr
c     do 200 qq=1,dimq
c     index(q)=qq
c     do 200 pp=1,dimp
c     index(p)=pp
c     b(index(1),index(2),index(3),index(4))=-a(pp,qq,rr,ss)
c200  continue
c
       if (p.eq.1) then
c     1***
       if (q.eq.2) then
c     12**
       if (r.eq.3) then
c     123* (4)
       do 211 ss=1,dims
       do 211 rr=1,dimr
       do 211 qq=1,dimq
       do 211 pp=1,dimp
       b(pp,qq,rr,ss)=-a(pp,qq,rr,ss)
 211    continue
       else
c     124* (3)
       do 212 ss=1,dims
       do 212 rr=1,dimr
       do 212 qq=1,dimq
       do 212 pp=1,dimp
       b(pp,qq,ss,rr)=-a(pp,qq,rr,ss)
 212    continue
       end if
       else if (q.eq.3) then
c     13**
       if (r.eq.2) then
c     132* (4)
       do 213 ss=1,dims
       do 213 rr=1,dimr
       do 213 qq=1,dimq
       do 213 pp=1,dimp
       b(pp,rr,qq,ss)=-a(pp,qq,rr,ss)
 213    continue
       else
c     134* (2)
       do 214 ss=1,dims
       do 214 rr=1,dimr
       do 214 qq=1,dimq
       do 214 pp=1,dimp
       b(pp,ss,qq,rr)=-a(pp,qq,rr,ss)
 214    continue
       end if
       else if (q.eq.4) then
c     14**
       if (r.eq.2) then
c     142* (3)
       do 215 ss=1,dims
       do 215 rr=1,dimr
       do 215 qq=1,dimq
       do 215 pp=1,dimp
       b(pp,rr,ss,qq)=-a(pp,qq,rr,ss)
 215    continue
       else
c     143* (2)
       do 216 ss=1,dims
       do 216 rr=1,dimr
       do 216 qq=1,dimq
       do 216 pp=1,dimp
       b(pp,ss,rr,qq)=-a(pp,qq,rr,ss)
 216    continue
       end if
       end if
c
       else if (p.eq.2) then
c     2***
       if (q.eq.1) then
c     21**
       if (r.eq.3) then
c     213* (4)
       do 221 ss=1,dims
       do 221 rr=1,dimr
       do 221 qq=1,dimq
       do 221 pp=1,dimp
       b(qq,pp,rr,ss)=-a(pp,qq,rr,ss)
 221    continue
       else
c     214* (3)
       do 222 ss=1,dims
       do 222 rr=1,dimr
       do 222 qq=1,dimq
       do 222 pp=1,dimp
       b(qq,pp,ss,rr)=-a(pp,qq,rr,ss)
 222    continue
       end if
       else if (q.eq.3) then
c     23**
       if (r.eq.1) then
c     231* (4)
       do 223 ss=1,dims
       do 223 rr=1,dimr
       do 223 qq=1,dimq
       do 223 pp=1,dimp
       b(rr,pp,qq,ss)=-a(pp,qq,rr,ss)
 223    continue
       else
c     234* (1)
       do 224 ss=1,dims
       do 224 rr=1,dimr
       do 224 qq=1,dimq
       do 224 pp=1,dimp
       b(ss,pp,qq,rr)=-a(pp,qq,rr,ss)
 224    continue
       end if
       else if (q.eq.4) then
c     24**
       if (r.eq.3) then
c     243* (1)
       do 225 ss=1,dims
       do 225 rr=1,dimr
       do 225 qq=1,dimq
       do 225 pp=1,dimp
       b(ss,pp,rr,qq)=-a(pp,qq,rr,ss)
 225    continue
       else
c     241* (3)
       do 226 ss=1,dims
       do 226 rr=1,dimr
       do 226 qq=1,dimq
       do 226 pp=1,dimp
       b(rr,pp,ss,qq)=-a(pp,qq,rr,ss)
 226    continue
       end if
       end if
c
       else if (p.eq.3) then
c     3***
       if (q.eq.1) then
c     31**
       if (r.eq.2) then
c     312* (4)
       do 231 ss=1,dims
       do 231 rr=1,dimr
       do 231 qq=1,dimq
       do 231 pp=1,dimp
       b(qq,rr,pp,ss)=-a(pp,qq,rr,ss)
 231    continue
       else
c     314* (2)
       do 232 ss=1,dims
       do 232 rr=1,dimr
       do 232 qq=1,dimq
       do 232 pp=1,dimp
       b(qq,ss,pp,rr)=-a(pp,qq,rr,ss)
 232    continue
       end if
       else if (q.eq.2) then
c     32**
       if (r.eq.1) then
c     321* (4)
       do 233 ss=1,dims
       do 233 rr=1,dimr
       do 233 qq=1,dimq
       do 233 pp=1,dimp
       b(rr,qq,pp,ss)=-a(pp,qq,rr,ss)
 233    continue
       else
c     324* (1)
       do 234 ss=1,dims
       do 234 rr=1,dimr
       do 234 qq=1,dimq
       do 234 pp=1,dimp
       b(ss,qq,pp,rr)=-a(pp,qq,rr,ss)
 234    continue
       end if
       else if (q.eq.4) then
c     34**
       if (r.eq.1) then
c     341* (2)
       do 235 ss=1,dims
       do 235 rr=1,dimr
       do 235 qq=1,dimq
       do 235 pp=1,dimp
       b(rr,ss,pp,qq)=-a(pp,qq,rr,ss)
 235    continue
       else
c     342* (1)
       do 236 ss=1,dims
       do 236 rr=1,dimr
       do 236 qq=1,dimq
       do 236 pp=1,dimp
       b(ss,rr,pp,qq)=-a(pp,qq,rr,ss)
 236    continue
       end if
       end if
c
       else if (p.eq.4) then
c     4***
       if (q.eq.1) then
c     41**
       if (r.eq.3) then
c     413* (2)
       do 241 ss=1,dims
       do 241 rr=1,dimr
       do 241 qq=1,dimq
       do 241 pp=1,dimp
       b(qq,ss,rr,pp)=-a(pp,qq,rr,ss)
 241    continue
       else
c     412* (3)
       do 242 ss=1,dims
       do 242 rr=1,dimr
       do 242 qq=1,dimq
       do 242 pp=1,dimp
       b(qq,rr,ss,pp)=-a(pp,qq,rr,ss)
 242    continue
       end if
       else if (q.eq.2) then
c     42**
       if (r.eq.1) then
c     421* (3)
       do 243 ss=1,dims
       do 243 rr=1,dimr
       do 243 qq=1,dimq
       do 243 pp=1,dimp
       b(rr,qq,ss,pp)=-a(pp,qq,rr,ss)
 243    continue
       else
c     423* (1)
       do 244 ss=1,dims
       do 244 rr=1,dimr
       do 244 qq=1,dimq
       do 244 pp=1,dimp
       b(ss,qq,rr,pp)=-a(pp,qq,rr,ss)
 244    continue
       end if
       else if (q.eq.3) then
c     43**
       if (r.eq.1) then
c     431* (2)
       do 245 ss=1,dims
       do 245 rr=1,dimr
       do 245 qq=1,dimq
       do 245 pp=1,dimp
       b(rr,ss,qq,pp)=-a(pp,qq,rr,ss)
 245    continue
       else
c     432* (1)
       do 246 ss=1,dims
       do 246 rr=1,dimr
       do 246 qq=1,dimq
       do 246 pp=1,dimp
       b(ss,rr,qq,pp)=-a(pp,qq,rr,ss)
 246    continue
       end if
       end if
c
       end if
c
       end if
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(s)
       end
c
c     --------------------------
c
       subroutine cct3_map31 (a,b,dimp,dimq,dimr,p,q,r,nfact)
c
c     maping A(p1,q1,r1) -> nfact*B(p2,q2,r2)
c
       real*8 a(*)
       real*8 b(*)
       integer dim(3)
       integer dimp,dimq,dimr,p,q,r,nfact
       dim(p)=dimp
       dim(q)=dimq
       dim(r)=dimr
       call cct3_map32 (a,b,dimp,dimq,dimr,dim(1),dim(2),dim(3),
     & p,q,r,nfact)
       return
       end
c
       subroutine cct3_map32 (a,b,dimp,dimq,dimr,dim1,dim2,dim3,
     & p,q,r,nfact)
c
       integer dimp,dimq,dimr,dim1,dim2,dim3,p,q,r,nfact
       real*8 a(1:dimp,1:dimq,1:dimr)
       real*8 b(1:dim1,1:dim2,1:dim3)
c     integer index(1:3)
c
       integer pp,qq,rr
c
       if (nfact.eq.1) then
c
c     fact = + 1
c
c     do 100 rr=1,dimr
c     index(r)=rr
c     do 100 qq=1,dimq
c     index(q)=qq
c     do 100 pp=1,dimp
c     index(p)=pp
c     b(index(1),index(2),index(3))=a(pp,qq,rr)
c100  continue
c
       if (p.eq.1) then
c     1**
       if (q.eq.2) then
c     12* (3)
       do 110 rr=1,dimr
       do 110 qq=1,dimq
       do 110 pp=1,dimp
       b(pp,qq,rr)=a(pp,qq,rr)
 110    continue
       else
c     13* (2)
       do 120 rr=1,dimr
       do 120 qq=1,dimq
       do 120 pp=1,dimp
       b(pp,rr,qq)=a(pp,qq,rr)
 120    continue
       end if
       else if (p.eq.2) then
c     2**
       if (q.eq.1) then
c     21* (3)
       do 130 rr=1,dimr
       do 130 qq=1,dimq
       do 130 pp=1,dimp
       b(qq,pp,rr)=a(pp,qq,rr)
 130    continue
       else
c     23* (1)
       do 140 rr=1,dimr
       do 140 qq=1,dimq
       do 140 pp=1,dimp
       b(rr,pp,qq)=a(pp,qq,rr)
 140    continue
       end if
       else if (p.eq.3) then
c     3**
       if (q.eq.1) then
c     31* (2)
       do 150 rr=1,dimr
       do 150 qq=1,dimq
       do 150 pp=1,dimp
       b(qq,rr,pp)=a(pp,qq,rr)
 150    continue
       else
c     32* (1)
       do 160 rr=1,dimr
       do 160 qq=1,dimq
       do 160 pp=1,dimp
       b(rr,qq,pp)=a(pp,qq,rr)
 160    continue
       end if
       end if
c
c
       else
c
c     factor = -1
c
c     do 200 rr=1,dimr
c     index(r)=rr
c     do 200 qq=1,dimq
c     index(q)=qq
c     do 200 pp=1,dimp
c     index(p)=pp
c     b(index(1),index(2),index(3))=-a(pp,qq,rr)
c     200        continue
c
       if (p.eq.1) then
c     1**
       if (q.eq.2) then
c     12* (3)
       do 210 rr=1,dimr
       do 210 qq=1,dimq
       do 210 pp=1,dimp
       b(pp,qq,rr)=-a(pp,qq,rr)
 210    continue
       else
c     13* (2)
       do 220 rr=1,dimr
       do 220 qq=1,dimq
       do 220 pp=1,dimp
       b(pp,rr,qq)=-a(pp,qq,rr)
 220    continue
       end if
       else if (p.eq.2) then
c     2**
       if (q.eq.1) then
c     21* (3)
       do 230 rr=1,dimr
       do 230 qq=1,dimq
       do 230 pp=1,dimp
       b(qq,pp,rr)=-a(pp,qq,rr)
 230    continue
       else
c     23* (1)
       do 240 rr=1,dimr
       do 240 qq=1,dimq
       do 240 pp=1,dimp
       b(rr,pp,qq)=-a(pp,qq,rr)
 240    continue
       end if
       else if (p.eq.3) then
c     3**
       if (q.eq.1) then
c     31* (2)
       do 250 rr=1,dimr
       do 250 qq=1,dimq
       do 250 pp=1,dimp
       b(qq,rr,pp)=-a(pp,qq,rr)
 250    continue
       else
c     32* (1)
       do 260 rr=1,dimr
       do 260 qq=1,dimq
       do 260 pp=1,dimp
       b(rr,qq,pp)=-a(pp,qq,rr)
 260    continue
       end if
       end if
c
       end if
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(r)
       end
c
c     --------------------------
c
       subroutine cct3_map21 (a,b,dimp,dimq,p,q,nfact)
c
c     maping A(p1,q1) -> nfact*B(p2,q2)
c
       real*8 a(*)
       real*8 b(*)
       integer dim(2)
       integer dimp,dimq,p,q,nfact
       dim(p)=dimp
       dim(q)=dimq
       call cct3_map22 (a,b,dimp,dimq,dim(1),dim(2),p,q,nfact)
       return
       end
c
       subroutine cct3_map22 (a,b,dimp,dimq,dim1,dim2,p,q,nfact)
c
       integer dimp,dimq,dim1,dim2,p,q,nfact
       real*8 a(1:dimp,1:dimq)
       real*8 b(1:dim1,1:dim2)
c       integer index(1:2)
c
       integer pp,qq
c
       if (nfact.eq.1) then
c
c     nfact = + 1
c
c     do 100 qq=1,dimq
c     index(q)=qq
c     do 100 pp=1,dimp
c     index(p)=pp
c     b(index(1),index(2))=a(pp,qq)
c100  continue
c
       if (p.eq.1) then
c     1* (2)
       do 110 qq=1,dimq
       do 110 pp=1,dimp
       b(pp,qq)=a(pp,qq)
 110    continue
       else
c     2* (1)
       do 120 qq=1,dimq
       do 120 pp=1,dimp
       b(qq,pp)=a(pp,qq)
 120    continue
       end if
c
c
       else
c
c     nfact = -1
c
c     do 200 qq=1,dimq
c     index(q)=qq
c     do 200 pp=1,dimp
c     index(p)=pp
c     b(index(1),index(2))=-a(pp,qq)
c200  continue
c
       if (p.eq.1) then
c     1* (2)
       do 210 qq=1,dimq
       do 210 pp=1,dimp
       b(pp,qq)=-a(pp,qq)
 210    continue
       else
c     2* (1)
       do 220 qq=1,dimq
       do 220 pp=1,dimp
       b(qq,pp)=-a(pp,qq)
 220    continue
       end if
c
       end if
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(q)
       end
c
c     --------------------------
c
       subroutine cct3_map11 (a,b,dimp,nfact)
c
c     maping A(p) -> nfact*B(p)
c
       integer dimp,nfact
       real*8 a(1:dimp)
       real*8 b(1:dimp)
c
       integer pp
c
       if (nfact.eq.1) then
c
       do 100 pp=1,dimp
       b(pp)=a(pp)
 100    continue
c
       else
c
       do 200 pp=1,dimp
       b(pp)=-a(pp)
 200    continue
c
       end if
c
       return
       end
c
